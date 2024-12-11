#!/usr/bin/env python3
import os
import sys
import configparser
import subprocess
import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
import multiprocessing
from functools import partial
import shutil

class BioinformaticsPipeline:
    def __init__(self, config_file):
        """初始化流程"""
        self.pipeline_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 使用configparser加载配置
        config = configparser.ConfigParser()
        config.read(config_file)
        self.config = config
        
        # 初始化状态文件路径
        self.status_file = None
        self.logger = None
        
    def _init_status_file(self, output_dir, sample):
        """初始化状态文件"""
        status_dir = os.path.join(output_dir, "pipeline_status")
        os.makedirs(status_dir, exist_ok=True)
        self.status_file = os.path.join(status_dir, f"{sample}.status.json")
        
        if not os.path.exists(self.status_file):
            status = {
                "QC": False,
                "Alignment": False,
                "SNV": False,
                "IS": False,
                "SV": False,
                "CNV": False,
                "RefCNV": False
            }
            with open(self.status_file, 'w') as f:
                json.dump(status, f, indent=2)
    
    def _get_status(self, step):
        """获取步骤执行状态"""
        if not self.status_file:
            return False
        
        with open(self.status_file) as f:
            status = json.load(f)
        return status.get(step, False)
    
    def _set_status(self, step, completed=True):
        """设置步骤执行状态"""
        if not self.status_file:
            return
            
        with open(self.status_file) as f:
            status = json.load(f)
        
        status[step] = completed
        
        with open(self.status_file, 'w') as f:
            json.dump(status, f, indent=2)
    
    def _should_run_step(self, step, force_all=False, skip_steps=None):
        """判断是否应该运行某个步骤"""
        # 检查是否在跳过列表中
        if skip_steps and step in skip_steps:
            self.logger.info(f"{step} 被指定跳过")
            return False
        
        if force_all:
            self.logger.info(f"{step} 强制执行")
            return True
        
        if self._get_status(step):
            self.logger.info(f"{step} 已完成")
            return False
        return True
    
    def _setup_logging(self, output_dir, sample):
        """设置日志记录"""
        log_dir = os.path.join(output_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        
        # 设置日志文件
        log_file = os.path.join(log_dir, f"{sample}.pipeline.log")
        
        # 配置日志格式
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Pipeline started for sample: {sample}")
    
    def _log_step_time(self, step, start=True):
        """记录步骤的开始或结束时间"""
        if start:
            self.logger.info(f"开始执行 {step}")
        else:
            self.logger.info(f"完成执行 {step}")
    
    def run_command(self, cmd, step_name):
        """执行命令并检查返回状态"""
        self._log_step_time(step_name, start=True)
        
        self.logger.info(f"执行命令: {cmd}")
        
        try:
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if process.returncode != 0:
                self.logger.error(f"{step_name} 执行失败:")
                self.logger.error(process.stderr)
                return False
                
            self._log_step_time(step_name, start=False)
            return True
            
        except Exception as e:
            self.logger.error(f"{step_name} 执行出错: {str(e)}")
            return False
    
    def get_adapter_sequences(self, platform=None):
        """获取指定平台的接头序列"""
        if platform is None:
            platform = self.config.get("adapters", "default_platform")
        
        platform = platform.lower()
        if platform == "illumina":
            return (self.config.get("adapters", "illumina_r1"),
                    self.config.get("adapters", "illumina_r2"))
        elif platform == "nextera":
            return (self.config.get("adapters", "nextera_r1"),
                    self.config.get("adapters", "nextera_r2"))
        elif platform == "bgi":
            return (self.config.get("adapters", "bgi_r1"),
                    self.config.get("adapters", "bgi_r2"))
        else:
            raise ValueError(f"Unsupported sequencing platform: {platform}")
    
    def run_qc(self, fq1, fq2, output_dir, sample, platform=None):
        """运行质控步骤"""
        if not self._should_run_step("QC"):
            return True
            
        # 获取对应平台的接头序列
        adapter1, adapter2 = self.get_adapter_sequences(platform)
        
        # 构建QC命令
        cmd = (f"bash {self.pipeline_dir}/QC.sh "
               f"{fq1} {fq2} "
               f"{adapter1} {adapter2} "
               f"{output_dir} {sample}")
        
        success = self.run_command(cmd, "QC")
        if success:
            self._set_status("QC")
        return success
    
    def run_align(self, fq1, fq2, output_dir, sample):
        """运行比对步骤"""
        ref_genome = self.config.get("paths", "genome")
        cmd = f"bash {self.pipeline_dir}/Align.sh {fq1} {fq2} {ref_genome} {output_dir} {sample}"
        success = self.run_command(cmd, "Alignment")
        if success:
            self._set_status("Alignment")
        return success
    
    def run_snv(self, bam, output_dir, sample):
        """运行SNV/Indel检测分析"""
        genome = self.config.get("paths", "genome")
        cmd = f"bash {self.pipeline_dir}/SNV.sh \
        {bam} \
        {genome} \
        {output_dir} \
        {sample}"
        
        return self.run_command(cmd, "SNV Analysis")
    
    def run_is_pipeline(self, bam, output_dir, sample):
        """运行IS分析流程"""
        ref_genome = self.config.get("paths", "genome")
        ltr_bed = self.config.get("paths", "ltr_bed")
        gene_exon_bed = self.config.get("paths", "gene_exon_bed")
        gtf_file = self.config.get("paths", "gtf_file")
        cmd = f"bash {self.pipeline_dir}/ISPipline.sh {bam} {ref_genome} {ltr_bed} {gene_exon_bed} {gtf_file} {output_dir} {sample}"
        return self.run_command(cmd, "IS Pipeline")
    
    def run_sv(self, bam, output_dir, sample):
        """运行结构变异分析"""
        ref_genome = self.config.get("paths", "genome")
        cmd = f"bash {self.pipeline_dir}/SV.sh {bam} {ref_genome} {output_dir} {sample}"
        return self.run_command(cmd, "SV Analysis")
    
    def run_cnv(self, bam, output_dir, sample):
        """运行拷贝数变异分析"""
        ref_genome = self.config.get("paths", "genome")
        genome_bed = self.config.get("paths", "genome_bed")
        window_size = self.config.get("parameters", "window_size")
        cmd = f"bash {self.pipeline_dir}/CNV.sh {bam} {ref_genome} {genome_bed} {window_size} {output_dir} {sample}"
        return self.run_command(cmd, "CNV Analysis")
    
    def run_reference_cnv(self, output_dir, sample, d4_file):
        """运行参考CNV分析"""
        ref_bed = self.config.get("paths", "ref_bed")
        target_bed = self.config.get("paths", "target_bed")
        cmd = f"bash {self.pipeline_dir}/reference.cnv.sh {d4_file} {ref_bed} {target_bed} {output_dir} {sample}"
        return self.run_command(cmd, "Reference CNV Analysis")
    
    def run_parallel_analysis(self, bam_file, args):
        """并行运行SNV、IS、SV、CNV和RefCNV分析"""
        try:
            # 创建进程池
            pool = multiprocessing.Pool(processes=5)  # 5个并行进程
            
            # 解析跳过步骤
            skip_steps = set(args.skip.split(',')) if args.skip else set()
            
            # 准备各个分析的参数
            analyses = []
            
            # SNV Analysis
            if self._should_run_step("SNV", args.force_all, skip_steps):
                snv_dir = os.path.join(args.output_dir, "SNV", args.sample)
                analyses.append(('SNV', partial(self.run_snv, bam_file, snv_dir, args.sample)))
            
            # IS Pipeline
            if self._should_run_step("IS", args.force_all, skip_steps):
                is_dir = os.path.join(args.output_dir, "IS", args.sample)
                analyses.append(('IS', partial(self.run_is_pipeline, bam_file, is_dir, args.sample)))
            
            # SV Analysis
            if self._should_run_step("SV", args.force_all, skip_steps):
                sv_dir = os.path.join(args.output_dir, "SV", args.sample)
                analyses.append(('SV', partial(self.run_sv, bam_file, sv_dir, args.sample)))
            
            # CNV Analysis
            if self._should_run_step("CNV", args.force_all, skip_steps):
                cnv_dir = os.path.join(args.output_dir, "CNV", args.sample)
                analyses.append(('CNV', partial(self.run_cnv, bam_file, cnv_dir, args.sample)))

            # RefCNV Analysis
            if self._should_run_step("RefCNV", args.force_all, skip_steps):
                ref_cnv_dir = os.path.join(args.output_dir, "RefCNV", args.sample)
                d4_file = os.path.join(args.output_dir, "Align", args.sample, f"{args.sample}.per-base.d4")
                analyses.append(('RefCNV', partial(self.run_reference_cnv, ref_cnv_dir, args.sample, d4_file)))
            
            # 并行执行分析
            results = []
            for name, func in analyses:
                self.logger.info(f"提交 {name} 分析任务")
                results.append((name, pool.apply_async(func)))
            
            # 等待所有任务完成并检查结果
            success = True
            for name, result in results:
                try:
                    if not result.get():  # 等待任务完成并获取结果
                        self.logger.error(f"{name} 分析失败")
                        success = False
                    else:
                        self._set_status(name)
                        self.logger.info(f"{name} 分析完成")
                except Exception as e:
                    self.logger.error(f"{name} 分析出错: {str(e)}")
                    success = False
            
            # 关闭进程池
            pool.close()
            pool.join()
            
            return success
            
        except Exception as e:
            self.logger.error(f"并行分析失败: {str(e)}")
            return False

    
    def get_ltr_result(self):
        """获取载体上的突变信息"""

        ltr_bed = self.config.get("paths", "ltr_bed")
        ltr_name = set()
        with open(ltr_bed, 'r') as f:
            for line in f:
                line_split = line.strip().split('\t')
                ltr_name.add(line_split[0])
        
        return ltr_name
    
    
    def _copy_results_to_result_dir(self, output_dir, sample):
        """将分析结果文件复制到Result目录"""
        # 创建结果目录
        result_dir = os.path.join(output_dir, "Result")
        os.makedirs(result_dir, exist_ok=True)    
        result_files = {
            # RefCNV结果
            os.path.join("RefCNV", f"{sample}.gene_average.txt"): 
                f"{sample}.target.copy.txt",
            
            # IS结果
            os.path.join("IS", sample, f"{sample}.insertion_break_point.anno.xls"): 
                f"{sample}.is.anno.xls" 
        }
        
        # 复制文件
        for src_path, dst_name in result_files.items():
            src_file = os.path.join(output_dir, src_path)
            dst_file = os.path.join(result_dir, dst_name)
            if os.path.exists(src_file):
                try:
                    shutil.copy2(src_file, dst_file)
                    self.logger.info(f"已复制结果文件: {dst_name}")
                except Exception as e:
                    self.logger.warning(f"复制文件失败 {src_file}: {str(e)}")

    def _summarize_vector_variants(self, output_dir, sample, ltr_names):
        """汇总载体上的变异信息"""
        result_dir = os.path.join(output_dir, "Result")
        os.makedirs(result_dir, exist_ok=True)
        
        # CNV结果文件
        cnv_variants_file = os.path.join(result_dir, f"{sample}.target.cnv.txt")
        with open(cnv_variants_file, 'w') as out:
            # 写入CNV表头
            header = ['Vector', 'Start', 'End', 'CNV_Type', 'Log2_Ratio']
            out.write('\t'.join(header) + '\n')
            
            # 处理CNV结果
            cnv_file = os.path.join(output_dir, "CNV", sample, f"{sample}_cnvkit_output", f"{sample}.segments.anno.txt")
            if os.path.exists(cnv_file):
                with open(cnv_file) as f:
                    next(f)  # 跳过表头
                    for line in f:
                        fields = line.strip().split('\t')
                        chr_name = fields[0]
                        if chr_name in ltr_names:
                            start = fields[1]
                            end = fields[2]
                            cnv_type = fields[3]
                            log2_ratio = fields[4]
                            out.write(f'{chr_name}\t{start}\t{end}\t{cnv_type}\t{log2_ratio}\n')
            else:
                self.logger.warning(f"CNV结果文件不存在: {cnv_file}")
        
        # SV结果文件
        sv_variants_file = os.path.join(result_dir, f"{sample}.target.sv.txt")
        with open(sv_variants_file, 'w') as out:
            # 写入SV表头
            header = ['Vector', 'Start', 'End', 'SV_Type', 'SV_ID']
            out.write('\t'.join(header) + '\n')
            
            # 处理SV结果
            sv_file = os.path.join(output_dir, "SV", sample, f"{sample}_sv.anno.xlsx")
            if os.path.exists(sv_file):
                with open(sv_file) as f:
                    next(f)  # 跳过表头
                    for line in f:
                        fields = line.strip().split('\t')
                        chr_name = fields[0]
                        if chr_name in ltr_names:
                            start = fields[1]
                            end = fields[2]
                            sv_type = fields[3]
                            sv_id = fields[4]
                            out.write(f'{chr_name}\t{start}\t{end}\t{sv_type}\t{sv_id}\n')
            else:
                self.logger.warning(f"SV结果文件不存在: {sv_file}")

        # SNV结果文件
        snv_variants_file = os.path.join(result_dir, f"{sample}.target.snv.txt")
        with open(snv_variants_file, 'w') as out:
            # 写入SNV表头
            header = ['Vector', 'Position', 'Variant_Type', 'Ref>Alt', 'Quality', 'Depth', 'Filter', 'AD', 'AF']
            out.write('\t'.join(header) + '\n')
            
            # 处理SNV结果
            snv_file = os.path.join(output_dir, "SNV", sample, f"{sample}.filtered.anno.xls")
            if os.path.exists(snv_file):
                with open(snv_file) as f:
                    next(f)  # 跳过表头
                    for line in f:
                        fields = line.strip().split('\t')
                        chr_name = fields[0]
                        if chr_name in ltr_names:
                            pos = fields[1]
                            var_type = fields[2]
                            ref_alt = fields[3]
                            qual = fields[4]
                            depth = fields[5]
                            filter_status = fields[6]
                            ad = fields[7]
                            af = fields[8]
                            out.write(f'{chr_name}\t{pos}\t{var_type}\t{ref_alt}\t{qual}\t'
                                    f'{depth}\t{filter_status}\t{ad}\t{af}\n')
            else:
                self.logger.warning(f"SNV结果文件不存在: {snv_file}")


    def run_pipeline(self, args):
        """运行完整的分析流程"""
        try:
            # 创建主输出目录
            os.makedirs(args.output_dir, exist_ok=True)
            
            # 设置日志记录
            self._setup_logging(args.output_dir, args.sample)
            
            # 解析跳过步骤
            skip_steps = set(args.skip.split(',')) if args.skip else set()
            
            # 记录运行参数
            self.logger.info("Pipeline parameters:")
            self.logger.info(f"  Sample: {args.sample}")
            self.logger.info(f"  Input FASTQ1: {args.fq1}")
            self.logger.info(f"  Input FASTQ2: {args.fq2}")
            self.logger.info(f"  Output directory: {args.output_dir}")
            self.logger.info(f"  Platform: {args.platform}")
            self.logger.info(f"  Force all: {args.force_all}")
            self.logger.info(f"  Skip steps: {args.skip if args.skip else 'None'}")
            
            # 初始化状态文件
            self._init_status_file(args.output_dir, args.sample)
            
            pipeline_start_time = datetime.now()
            self.logger.info(f"Pipeline execution started at: {pipeline_start_time}")
            
            # 1. QC
            qc_dir = os.path.join(args.output_dir, "QC", args.sample)
            if self._should_run_step("QC", args.force_all, skip_steps):
                if not self.run_qc(args.fq1, args.fq2, qc_dir, args.sample, args.platform):
                    self.logger.error("QC步骤失败，流程终止")
                    return False
            
            # 2. Alignment
            align_dir = os.path.join(args.output_dir, "Align", args.sample)
            if self._should_run_step("Alignment", args.force_all, skip_steps):
                fq1 = os.path.join(qc_dir, '{}.trim.R1.fastq.gz'.format(args.sample))
                fq2 = os.path.join(qc_dir, '{}.trim.R2.fastq.gz'.format(args.sample))
                if not self.run_align(fq1, fq2, align_dir, args.sample):
                    self.logger.error("Alignment步骤失败，流程终止")
                    return False
            
            bam_file = os.path.join(align_dir, '{}.deduped.bam'.format(args.sample))
            
            # 3. 并行运行所有后续分析
            if not self.run_parallel_analysis(bam_file, args):
                self.logger.error("并行分析步骤失败")
                return False

            # 4. 复制结果文件到Result目录
            self._copy_results_to_result_dir(args.output_dir, args.sample)
            
            # 5. 获取载体���称列表
            ltr_names = self.get_ltr_result()
            
            # 6. 汇总载体上的变异
            self._summarize_vector_variants(args.output_dir, args.sample, ltr_names)

            pipeline_end_time = datetime.now()
            duration = pipeline_end_time - pipeline_start_time
            self.logger.info(f"Pipeline execution completed at: {pipeline_end_time}")
            self.logger.info(f"Total execution time: {duration}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline execution failed: {str(e)}")
            return False
    


def main():
    pipeline_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    config_file = os.path.join(pipeline_dir, "config.ini")
    parser = argparse.ArgumentParser(description="Bioinformatics Analysis Pipeline")
    parser.add_argument("--fq1", required=True, help="Input fastq1 file")
    parser.add_argument("--fq2", required=True, help="Input fastq2 file")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--platform", default='illumina', choices=["illumina", "nextera", "bgi"],
                      help="Sequencing platform (default: illumina)")
    parser.add_argument("--force-all", action="store_true", 
                      help="强制执行所有步骤，忽略已完成状态")
    parser.add_argument("--skip", help="跳过指定步骤(逗号分隔): QC,Alignment,SNV,IS,SV,CNV,RefCNV")
    args = parser.parse_args()
    # 初始化并运行流程
    pipeline = BioinformaticsPipeline(config_file)
    if not pipeline.run_pipeline(args):
        sys.exit(1)

if __name__ == "__main__":
    main() 