using Pipe
using DataFrames

abstract PlinkOps
abstract PlinkInput <: PlinkOps

type PlinkBinaryFile <: PlinkInput
	resolved::Bool
	prev::PlinkOps
	next::Array{PlinkOps,1}
	bed::ASCIIString
	bim::ASCIIString
	fam::ASCIIString
	PlinkBinaryFile(bed,bim,fam) = (
		x = new();
		x.resolved=false;
		x.prev=x;
		x.next=[];
		x.bed=bed;
		x.bim=bim;
		x.fam=fam;
		x)

PlinkBinaryFile(p,bed,bim,fam) = (
		x = new();
		x.resolved=false;
		x.prev=p;
		x.next=[];
		x.bed=bed;
		x.bim=bim;
		x.fam=fam;
		x)
end

type PlinkFile <: PlinkInput
	resolved::Bool
	prev::PlinkOps
	next::Array{PlinkOps,1}
	ped::ASCIIString
	map::ASCIIString
	PlinkFile(ped,map) = (
		x = new();
		x.resolved=false;
		x.prev=x;
		x.next=[];
		x.ped=ped;
		x.map=map;
		x)
end

type PlinkVCF <: PlinkInput
	resolved::Bool
	prev::PlinkOps
	next::Array{PlinkOps,1}
	vcf::ASCIIString
end

type PlinkOption <: PlinkOps 
	resolved::Bool
	prev::PlinkOps
	next::Array{PlinkOps,1}
	op::Symbol
	args::Array{ASCIIString,1}
	files::Dict{Symbol,ASCIIString}
end

# Binary (or greater) PlinkOption with output files
function PlinkOption(prev::PlinkOps,op::Symbol,args::Array{ASCIIString,1},files::Dict{Symbol,ASCIIString}) 
	opt=PlinkOption(false,prev,Array{PlinkOps,1}(),op,args,files)
	push!(opt.prev.next,opt)
	opt
end	

# Binary (or greater) PlinkOption without output files
PlinkOption(prev::PlinkOps,op::Symbol,args::Array{ASCIIString,1}) = PlinkOption(prev,op,args,Dict{Symbol,ASCIIString}())

# Unary PlinkOption with no extra tokens or output files
PlinkOption(prev::PlinkOps,op::Symbol) = PlinkOption(prev,op,Array{ASCIIString,1}(),Dict{Symbol,ASCIIString}())

# Unary PlinkOption with file output with output files
PlinkOption(prev::PlinkOps,op::Symbol,files::Dict{Symbol,ASCIIString}) = PlinkOption(prev,op,Array{ASCIIString,1}(),files)

# Checks first for a binary with given plinkstem
# If not found, checks for a non-binary plink file
function plink(plinkstem::ASCIIString)
	bfile = map(x -> string(plinkstem,".",x), ["bed","bim","fam"])
	file = map(x -> string(plinkstem,".",x), ["ped","map"])
	vcf = string(plinkstem, ".vcf")

	if all(isfile,bfile)
		PlinkBinaryFile(bfile...)
	elseif all(isfile,file)
		PlinkFile(file...)
	elseif isfile(vcf)
		PlinkVCF(vcf)
	else
		error(string("No suitable plink files found for ", plinkstem))
	end
end

#plink(plinkstem::ASCIIString) = plink(ASCIIString(plinkstem))

sym_to_arg(s::Symbol) = string("--", replace(string(s),"_","-"))

function plink_arguments(p::PlinkInput)
	fields = setdiff(fieldnames(p),[:prev,:next,:resolved])
	opt = map(sym_to_arg,fields)
	arg = map(i->getfield(p,i), fields)
	paired = map(x->[x[1],x[2]],zip(opt,arg))
	foldl(append!,[],paired)
end

plink_arguments(p::PlinkOption) = vcat(sym_to_arg(p.op),p.args)

collect_until_input(p::PlinkInput) = [ p ] 
collect_until_input(p::PlinkOps) = [ p ; collect_until_input(p.prev) ]

function split_operations(p::PlinkOps)
	# A recode / make_bed operation always result
	# in a PlinkInput object, but it is not part of
	# the operations, but rather the result
	ptr = (typeof(p) <: PlinkInput)	? p.prev : p

	commands = Array{PlinkOps,1}()

	while ptr != ptr.prev 
		opseq = collect_until_input(ptr)
		push!(commands, ptr)
		ptr = last(opseq).prev
	end

	return commands
end

build_plink_commands(p::PlinkInput) = `plink $(plink_arguments(p))`
build_plink_commands(p::PlinkOps) = `$(build_plink_commands(p.prev))  $(plink_arguments(p))`

plinkcmds(p::PlinkOps) = map(x->build_plink_commands(x),split_operations(p))

#build_command(pieces::Array{ASCIIString,1}) = foldl((x,y)->`$x $y`,``,pieces)

function get_file_name(p::PlinkOps, ext::ASCIIString)
	# Basecase: If there is no specification of out-name, so just assume "plink"
	if (typeof(p) <: PlinkInput)
		string("plink",ext)
	elseif (typeof(p)==PlinkOption && p.op == :out)
		string(p.args[1],ext) 
	else
		get_file_name(p.prev,ext)
	end
end


function out(p::PlinkOps,dirpath::ASCIIString)
	# FIXME: Dirpath could be checked to be valid
	# FIXME: Check that a previous out was not allready specified 
	PlinkOption(p,:out,[dirpath])
end

# Create "anonymous" plink option 
function opt(p::PlinkOps,opt::Cmd)
	m=match(r"`--(\S+)\s?(.*)`",string(opt))
	println(m)
	if m[1]==nothing
		error(string("Invalid unchecked plink option: ", opt))
	elseif m[2]==""
		PlinkOption(p,Symbol(m[1]))
	else
		PlinkOption(p,Symbol(m[1]),map(ASCIIString,split(reverse(chomp(reverse(m[2])))," ")))
	end
end

extensions_dict(p,names::Array{Symbol,1}) = [ x => ASCIIString(get_file_name(p,string(".",x))) for x in names ]

##################################################################################
## Input filtering
##################################################################################

macro file_filter_function(name,namesym)
	eval(:(($name)(p::PlinkOps,f::ASCIIString ; range=false) = 
		(isfile(f) ? 
			PlinkOption(p,$namesym,range?["range",f]:[f])
			:
			error(string("No such file '", f, "' for plink option"))
		)
	))
end

## ID lists:
""" 
`keep(p::PlinkOption,file)` accepts a space/tab-delimited text `file` with family 
IDs in the first column and within-family IDs in the second column, 
and removes all unlisted samples from the current analysis. 
`remove` does the same for all listed samples.
Similarly `keep(p::PlinkOption,listed::Array{ASCIIString,2})` accepts a 2d array 
and autogenerates `file` for plink.
"""
@file_filter_function(keep,:keep)

""" 
`remove` accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, 
and removes all listed samples from the current analysis. `keep` does the same for all unlisted samples.
"""
@file_filter_function(exclude,:exclude)


""" 
`extract` normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.
"""
@file_filter_function(extract,:extract)

"""`from` excludes all variants on different chromosomes than the named variant, as well as those with smaller base-pair position values."""
from(p::PlinkOps,variant) = PlinkOption(p,:from, [ variant ])

"""
--from excludes all variants on different chromosomes than the named variant, as well as those with smaller base-pair position values. --to is similar, excluding variants with larger position values instead. If they are used together but the --from variant is after the --to variant, they are automatically swapped.
"""
to(p::PlinkOps,variant::ASCIIString) = PlinkOption(p,:to, [ variant ])

snp(p::PlinkOps,variant::ASCIIString) = PlinkOption(p,:snp, [ variant ])

window(p::PlinkOps,windows_size_kb::Int32) = PlinkOption(p,:snp, [ string(windows_size_kb) ])

exclude_snp(p::PlinkOps,variant::ASCIIString) = PlinkOption(p,:exclude_snp, [ variant ])

""" 
geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed
"""
geno(p::PlinkOps,threshold::Float64) = PlinkOption(p,:geno, [ string(threshold) ])

""" 
mind filters out all individuals with missing call rates exceeding the provided value (default 0.1) to be removed
"""
mind(p::PlinkOps,threshold::Float64) = PlinkOption(p,:geno, [ string(threshold) ])


maf(p::PlinkOps,threshold::Float64) = PlinkOption(p,:maf, [ string(threshold) ]) 
max_maf(p::PlinkOps,threshold::Float64) = PlinkOption(p,:max_maf, [ string(threshold) ]) 

mac(p::PlinkOps,threshold::Int32) = PlinkOption(p,:mac, [ string(threshold) ]) 
max_mac(p::PlinkOps,threshold::Int32) = PlinkOption(p,:max_mac, [ string(threshold) ]) 


"""filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold."""
hwe(p::PlinkOps,threshold::Float64) = PlinkOption(p,:hwe, [ string(threshold) ]) 

""" filters out all samples with missing phenotypes """
prune(p::PlinkOps) = PlinkOption(p,:prune) 

filter_cases(p::PlinkOps) = PlinkOption(p,:filter_cases) 
filter_controls(p::PlinkOps) = PlinkOption(p,:filter_controls) 
filter_males(p::PlinkOps) = PlinkOption(p,:filter_males) 
filter_females(p::PlinkOps) = PlinkOption(p,:filter_females) 
filter_founders(p::PlinkOps) = PlinkOption(p,:filter_founders) 
filter_nonfounders(p::PlinkOps) = PlinkOption(p,:filter_nonfounders) 

##################################################################################
## Basic statistics 
##################################################################################

freq(p::PlinkOps) = PlinkOption(p,:freq, Dict(:frq => get_file_name(p,".frq")))

hardy(p::PlinkOps) = PlinkOption(p,:hardy, Dict(:hwe => get_file_name(p,".hwe")))

ibc(p::PlinkOps) = PlinkOption(p,:ibc, Dict(:hwe => get_file_name(p,".ibc")))

check_sex(p::PlinkOps) = PlinkOption(p,:check_sex, Dict(:hwe => get_file_name(p,".sexcheck")))

#missing(p::PlinkOps,threshold::Float64) = PlinkOption(p,:missing, [ string(threshold) ]) 

##################################################################################
## Data management 
##################################################################################

function make_bed(p::PlinkOps)
	po = PlinkOption(p,:make_bed, extensions_dict(p,[:bed,:bim,:fam]))
	PlinkBinaryFile(po,po.files[:bed],po.files[:bim],po.files[:fam])
end

##################################################################################
## Reading and writing of plink file formats 
##################################################################################

# Basic types for representing markers
# Compared to the binary plinkformat, this is a quite wasteful representation.
# So it will take up more memory and it does when stored on disk 
# We may consider different options, e.g., runlength coding

type BiallelicVariantInfo
	chr::ASCIIString
	name::ASCIIString
	genetic_distance::Int32
	pos::Int32
	a1::ASCIIString # major allele
	a2::ASCIIString # minor allele

	# These are necessary to determine minor allele ordering
	a1_count::Int32
	a2_count::Int32
end

type BiallelicVariant
	info::BiallelicVariantInfo
	genotype::Int8
end

function Base.show(io::IO,x::BiallelicVariant) 
	if (x.genotype==0) 
		print(io, string(x.info.a1,"/",x.info.a1))
	elseif (x.genotype==1) 
		print(io, string(x.info.a1,"/",x.info.a2))
	else
		print(io, string(x.info.a1,"/",x.info.a2))
	end
end

function read_map_file(mapfile)
	markers_info = Array{BiallelicVariantInfo,1}()
	open(mapfile) do mapf
		for line in eachline(mapf)
			chr,snpname,cm,pos = split(chomp(line),[' ','\t'])
			push!(markers_info,BiallelicVariantInfo(chr,snpname,parse(Int32,cm),parse(Int32,pos),"","",0,0))
		end
	end
	return markers_info
end

function dataframe(p::PlinkFile) 
	markers_info = read_map_file(p.map)

	markers = Array{DataArray{BiallelicVariant,1},1}()

	family_id = data(Array{ASCIIString,1}()) 
	individual_id = data(Array{ASCIIString,1}()) 
	paternal_id = data(Array{ASCIIString,1}())
	maternal_id = data(Array{ASCIIString,1}())
	sex = data(Array{Int8,1}())
	phenotype = data(Array{Int8,1}())

	open(p.ped) do ped
		for line in eachline(ped) 
			fields = split(chomp(line),[' ','\t'])
			push!(family_id,shift!(fields))
			push!(individual_id,shift!(fields))
			push!(paternal_id,shift!(fields))
			push!(maternal_id,shift!(fields))
			push!(sex, parse(Int8,shift!(fields)))
			push!(phenotype, parse(Int8,shift!(fields)))
			for i in 1:length(markers_info)
				a1=shift!(fields)
				a2=shift!(fields)

				if markers_info[i].a1 == "" && a1 != "0"
					markers_info[i].a1 = a1	
				end

				if markers_info[i].a2 == "" && a2 != "0"
					markers_info[i].a2 = a2
				end

				genotype = if a1 == markers_info[i].a1 && a2 == markers_info[i].a1 		# Homozygous wildtype
					BiallelicVariant(markers_info[i],0)
				elseif a1 == markers_info[i].a2 && a2 == markers_info[i].a2		# Homozygous variant 
					BiallelicVariant(markers_info[i],2)
				elseif (a1!=a2) && (a1 == markers_info[i].a1) && (a2 == markers_info[i].a2)
					BiallelicVariant(markers_info[i],1)
				elseif (a1!=a2) && (a1 == markers_info[i].a2) && (a2 == markers_info[i].a1)
					BiallelicVariant(markers_info[i],1)
				else # in case of invalid or missing variant
					BiallelicVariant(markers_info[i],-9)
				end

				if length(markers) < i
					push!(markers,DataArray{BiallelicVariant,1}[])
				end

				push!(markers[i],genotype) 
				if genotype.genotype == -9
					markers[i][end] = NA
				end
			end
		end
	end

	df = DataFrame()
	df[:FID] = family_id
	df[:IID] = individual_id
	df[:PaternalID] = paternal_id 
	df[df[:PaternalID].==-9] = NA
	df[:MaternalID] = maternal_id 
	df[df[:MaternalID].==-9] = NA
	df[:Sex] = sex
	df[df[:Sex].==-9] = NA
	df[:Pheno] = phenotype
	df[df[:Pheno].==-9] = NA

	for i in 1:length(markers)
		df[Symbol(markers_info[i].name)] = markers[i] 
	end

	return df
end

dataframe(p::PlinkBinaryFile) = error("Not implemented yet")
dataframe(p::PlinkVCF) = error("Not implemented yet")

function dataframe(p::PlinkOps)
	for cmd in reverse(plinkcmds(p))
		run(cmd)
	end

	println(p.files)

	if (length(p.files) == 1)
		readtable(p.files[first(keys(p.files))],header=true,separator=' ')
	else
		error(string("Ambigous set of files-type: ", keys(p.files)))
	end
end

#function dataframe(p::PlinkOps,fileid::Symbol)
#	run(cmd(p))
#	if(p.files
#end
