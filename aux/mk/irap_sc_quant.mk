#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================


#************************************
# scRNA-seq - digital gene expression
#************************************
# UMI quantification can be made at gene or transcript level

# what tag shoud be used to get gene/transcript id
get_bam_tag=$(if $(filter $(1),gene),GX,tx)

# salmon+umi
# salmon - generate BAM
# UMI - count transcripts
#ifeq (salmon_umi,$(quant_method))
#salmon_quant_params+=--writeMappings=<outfile> 
#endif


## Kallisto+UMI
ifeq (kallisto_umi,$(quant_method))

# use kallisto to align

define make-kallisto-quant-rule=

$(call lib2quant_folder,$(1))$(1)/$(1).bam: $(call libname2fastq,$(1)) $(kallisto_index)
	(mkdir -p $$(@D) && $$(call run_kallisto_aln,$$(@D),$(call libname2fastq,$(1)),$($(1)_rs),$(call is_pe_lib,$(1))) && mv $$(@D)/abundance.tsv $$@)


$(call lib2quant_folder,$(1))$(2).transcripts.raw.kallisto.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,4 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@

$(call lib2quant_folder,$(1))$(2).transcripts.tpm.kallisto.kallisto.tsv: $(call lib2quant_folder,$(1))$(1)/$(1).abundance.tsv
	cut -f 1,5 $$< |tail -n +2 > $$@.tmp && mv $$@.tmp $$@


$(call lib2quant_folder,$(1))$(2).genes.raw.kallisto.tsv: $(call lib2quant_folder,$(1))$(2).transcripts.raw.kallisto.tsv $(mapTrans2gene)
	libTSVAggrTransByGene -i $$< -m $(mapTrans2gene) -o $$@.tmp && mv  $$@.tmp $$@


endef


ifneq ($(deps_check),nocheck)
# quantification
$(foreach l,$(se),$(eval $(call make-kallisto-quant-rule,$(l),$(l).se)))
$(foreach l,$(pe),$(eval $(call make-kallisto-quant-rule,$(l),$(l).pe)))
endif

endif



#************
# umis
#************

##$($(1)_known_umi_file)
ifeq (umis,$(quant_method))
# $1 - lib
# $2 - bam file prefix (includes .se|.pe)
# $3 - gene|transcript
# $4 - column in the reference file (mapTrans2gene) with all genes/transcripts
define make-umis-rule=
# 
$(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).mtx.gz: $(call lib2bam_folder,$(1))$(2).hits.bam  $$(mapTrans2gene) 
	mkdir -p $$(@D) && \
	time umis tagcount $$(umis_params) --sparse --parse_tags --gene_tags $$< $$@.tmp   && \
	$(if $($(1)_sample_name),sed "s/$$$$/$($(1)_sample_name)/" $$@.tmp.colnames,cat $$@.tmp.colnames) | gzip -f -c -  > $$(subst .gz,,$$@)_cols.gz && \
	gzip -f -c $$@.tmp.rownames > $$(subst .gz,,$$@)_rows.gz && \
	gzip -f -c $$@.tmp > $$@ && \
	rm -f $$@.tmp.* || ( rm -f $$@* && exit 1)
endef
#	time umis tagcount $$(umis_params) --sparse --parse_tags --gene_tags $$< /dev/stdout | tr "," "\t" > $$@.tmp   && \
#	add_missing_features --sort --tsv $$@.tmp --all_feat $$(mapTrans2gene) --all_feat_col $(4) --out $$@.tmp2 &&\

ifneq ($(deps_check),nocheck)
# gene level quantification
$(foreach l,$(se),$(eval $(call make-umis-rule,$(l),$(l).se,gene,1)))	
$(foreach l,$(pe),$(eval $(call make-umis-rule,$(l),$(l).pe,gene,1)))
endif

ifeq ($(transcript_expr),y)
$(call p_error,Unable to get transcript level quantification while using "umis". Please select another value for the transcript_quant option.)
endif

# Generate a single matrix
$(quant_toplevel_folder)/genes.raw.$(quant_method).mtx.gz: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).mtx.gz) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).mtx.gz) 
	( $(call pass_args_stdin,irap_merge_mtx,$@,-o $@ --all_feat $(mapTrans2gene) --all_feat_col 1 --in "$^") ) || ( rm -f $@ && exit 1)


endif
## umis


#*****************
# simple_umi_count
#*****************

ifeq (umi_count,$(quant_method))

#WAVE3_p_TARGETS+=$(subst .hits.bam,.hits.bytag_$(CELL_TAG).bam,$(WAVE2_TARGETS))


# $1 - lib
# $2 - bam file prefix (includes .se|.pe)
# $3 - gene|transcript
# $4 - column in the reference file (mapTrans2gene) with all genes/transcripts
define make-iumi-count-rule=
$(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).mtx.gz: $(call lib2bam_folder,$(1))$(2).hits.bytag_$(CELL_TAG).bam  $$(mapTrans2gene) $($(1)_known_umi_file) $($(1)_known_cells_file)
	time bam_umi_count --sorted_by_cell --cell_tag $(CELL_TAG) --bam $$< --ucounts $(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).mtx --rcounts $(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).reads.mtx $$(bam_umi_count_params) $(if $($(1)_known_umi_file),--known_umi $($(1)_known_umi_file),) $(if $($(1)_known_cells_file),--known_cells $($(1)_known_cells_file),) --tag $(call get_bam_tag,$(3))  $(if $($(1)_sample_name),--cell_suffix "-$($(1)_sample_name)",) --max_cells $(sc_max_cells) --max_feat $(sc_max_features) --feat_cell $(sc_feat_cell) && gzip -f $(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).mtx_{rows,cols} && gzip -f $(call lib2quant_folder,$(1))$(2).$(3)s.raw.$(quant_method).mtx || ( rm -f $$@ && exit 1)
endef

# irap_sc conf=conf/sc/10x_v2_pbmc.conf stage3  mapper=bowtie2 quant_method=irap_umi_count se=SE1

#irap_sc conf=conf/sc/10x_v2_pbmc.conf stage2  mapper=bowtie2 stage3 quant_method=irap_umi_count

ifneq ($(deps_check),nocheck)
# gene level quantification
$(foreach l,$(se),$(eval $(call make-iumi-count-rule,$(l),$(l).se,gene,1)))	
$(foreach l,$(pe),$(eval $(call make-iumi-count-rule,$(l),$(l).pe,gene,1)))
endif

ifeq ($(transcript_expr),y)
ifneq ($(deps_check),nocheck)
$(foreach l,$(se),$(eval $(call make-iumi-count-rule,$(l),$(l).se,transcript,2)))	
$(foreach l,$(pe),$(eval $(call make-iumi-count-rule,$(l),$(l).pe,transcript,2)))
endif
## 
$(quant_toplevel_folder)/transcripts.raw.$(quant_method).mtx.gz: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.transcripts.raw.$(quant_method).mtx.gz) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.transcripts.raw.$(quant_method).mtx.gz)
	( $(call pass_args_stdin,irap_merge_mtx,$@, -o $@ --all_feat $(mapTrans2gene) --all_feat_col 2 --in "$^") ) || ( rm -f $@ && exit 1)
endif

# Generate a single matrix
$(quant_toplevel_folder)/genes.raw.$(quant_method).mtx.gz: $(foreach p,$(pe),$(call lib2quant_folder,$(p))$(p).pe.genes.raw.$(quant_method).mtx.gz) $(foreach s,$(se), $(call lib2quant_folder,$(s))$(s).se.genes.raw.$(quant_method).mtx.gz) 
	( $(call pass_args_stdin,irap_merge_mtx,$@, -o $@ --all_feat_col 1 --all_feat $(mapTrans2gene) --in "$^") ) || ( rm -f $@ && exit 1)


endif
## umi_count

