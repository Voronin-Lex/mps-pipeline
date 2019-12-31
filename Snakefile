
def files():
    list = []
    for entry in os.listdir('./input'):
        list.append("{}".format(entry.split('.')[0]))
    return list

def create_directories():    
    os.makedirs("output/reports", exist_ok=True)
    os.makedirs("output/freebayes", exist_ok=True)
    os.makedirs("output/haplotypecaller", exist_ok=True)
    os.makedirs("output/samtools", exist_ok=True)
    os.makedirs("output/intersected", exist_ok=True)
    os.makedirs("output/individual", exist_ok=True)
    os.makedirs("output/common", exist_ok=True)

def gen():
    list = []
    for entry in os.listdir('./input'):
        list.append("output/reports/{}".format(entry.split('.')[0]))
    return list
    
create_directories()

samples = files()

terminate = gen()

rule terminate:
  input: terminate

rule countall:
  input:
    "output/individual/{sample}",
    "output/common/{sample}"
  output:
    "output/reports/{sample}"
  run:      
      with open(input[0], "r") as src1:
        vc_2 = src1.readline().split()[-2]
        vc_3 = src1.readline().split()[-2]
        vc_1 = src1.readline().split()[-2]
      with open(input[1], "r") as src2:
        vc_12 = src2.readline().split()[-2]
        vc_13 = src2.readline().split()[-2]
        vc_23 = src2.readline().split()[-2]
      with open(output[0], "w") as out:
      	out.write("{}\t{}\t{}\n".format(vc_1, vc_12, vc_13))
        out.write("{}\t{}\t{}\n".format(vc_2, "None", vc_23))
        out.write("{}\t{}\t{}".format(vc_3, "None", "None"))

rule countcommon:
  input:
    "output/intersected/{sample}_vc_1_vc_2.vcf.gz",
    "output/intersected/{sample}_vc_1_vc_3.vcf.gz",
    "output/intersected/{sample}_vc_2_vc_3.vcf.gz"
  output:
    "output/common/{sample}"
  shell:
    "vcftools --gzvcf {input[0]} 2>&1 >/dev/null | grep 'a possible' >> {output}; "
    "vcftools --gzvcf {input[1]} 2>&1 >/dev/null | grep 'a possible' >> {output}; "
    "vcftools --gzvcf {input[2]} 2>&1 >/dev/null | grep 'a possible' >> {output}; "


rule countvars:
  input:    
    "output/haplotypecaller/{sample}.vcf.gz",
    "output/samtools/{sample}.vcf.gz",
    "output/freebayes/{sample}.vcf.gz"
  output:
    "output/individual/{sample}"
  shell:
    "vcftools --gzvcf {input[0]} 2>&1 >/dev/null | grep 'a possible' >> {output}; "
    "vcftools --gzvcf {input[1]} 2>&1 >/dev/null | grep 'a possible' >> {output}; "
    "vcftools --gzvcf {input[2]} 2>&1 >/dev/null | grep 'a possible' >> {output}; " 


rule intersectvars:
  input:
    "output/freebayes/{sample}.vcf.gz",
    "output/haplotypecaller/{sample}.vcf.gz",
    "output/samtools/{sample}.vcf.gz",
    "output/freebayes/{sample}.vcf.gz.tbi",
    "output/haplotypecaller/{sample}.vcf.gz.tbi",
    "output/samtools/{sample}.vcf.gz.tbi"
  output:
    "output/intersected/{sample}_vc_1_vc_2.vcf.gz",
    "output/intersected/{sample}_vc_1_vc_3.vcf.gz",
    "output/intersected/{sample}_vc_2_vc_3.vcf.gz"
  shell:
    "vcf-isec -f -n +2 {input[0]} {input[1]} | bgzip -c > {output[0]}; "
    "vcf-isec -f -n +2 {input[0]} {input[2]} | bgzip -c > {output[1]}; "
    "vcf-isec -f -n +2 {input[1]} {input[2]} | bgzip -c > {output[2]}"


rule tabix:
  input:
    "output/haplotypecaller/{sample}.vcf.gz",
    "output/samtools/{sample}.vcf.gz",
    "output/freebayes/{sample}.vcf.gz"
  output:
    "output/haplotypecaller/{sample}.vcf.gz.tbi",
    "output/samtools/{sample}.vcf.gz.tbi", 
    "output/freebayes/{sample}.vcf.gz.tbi" 
  shell:
    "tabix -p vcf {input[0]}; "
    "tabix -p vcf {input[1]}; "
    "tabix -p vcf {input[2]}"

rule indexbgzip:
  input: 
    "output/haplotypecaller/{sample}.vcf",
    "output/samtools/{sample}.vcf",
    "output/freebayes/{sample}.vcf"
  output: 
    "output/haplotypecaller/{sample}.vcf.gz",
    "output/samtools/{sample}.vcf.gz", 
    "output/freebayes/{sample}.vcf.gz", 
  shell:
    "bgzip {input[0]}; "
    "bgzip {input[1]}; "
    "bgzip {input[2]}"

rule haplotypecaller:
  input: 
    "data/22.fa",
    "data/22.fa.fai",
    "data/22.dict",	
    expand("input/{sample}.bam", sample = samples),
    expand("input/{sample}.bam.bai", sample = samples)
  output:
    "output/haplotypecaller/{sample}.vcf"   
  shell: 
    "java -jar $GATK -R {input[0]} -T HaplotypeCaller -I {input[3]} -o {output}"

rule bamidx:
  input:
    "input/{sample}.bam"
  output:
    "input/{sample}.bam.bai"
  shell:
    "samtools index {input}" 


rule refdict:
  input:
    "data/22.fa"
  output:
    "data/22.dict"
  shell:
    "java -jar $PICARD CreateSequenceDictionary R={input} O={output}"

rule refidx:
  input:
    "data/22.fa"
  output:
    "data/22.fa.fai"
  shell:
    "samtools faidx {input}"

rule samtools:
  input:
    "data/22.fa",
    expand("input/{sample}.bam", sample = samples)
  output:
    "output/samtools/{sample}.vcf"
  shell:  
    "samtools mpileup -uf {input[0]} {input[1]} | bcftools view -vcg - > {output}"


rule freebayes:
  input:
    "data/22.fa",
    expand("input/{sample}.bam", sample = samples)
  output:
    "output/freebayes/{sample}.vcf"
  shell:
    "freebayes -f {input[0]} {input[1]} > {output}"



