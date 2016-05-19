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
	elseif (x.genotype==2)
		print(io, string(x.info.a2,"/",x.info.a2))
	else
		"??"
	end
end

each_trimmed_nonempty_line(io::IOStream) = filter(l->l!="", map(l->chomp(lstrip(l)), eachline(io)))
getfields(line) = filter(x->x!="",split(chomp(lstrip(line)),[' ','\t']))

function map_readfile(mapfile)
	println(readall(mapfile))
	markers_info = Array{BiallelicVariantInfo,1}()
	open(mapfile) do mapf
		for line in each_trimmed_nonempty_line(mapf)
			chr,snpname,cm,pos = getfields(line)
			push!(markers_info,BiallelicVariantInfo(chr,snpname,parse(Int32,cm),parse(Int32,pos),"","",0,0))
		end
	end
	return markers_info
end

function bim_readfile(bimfile)
	markers_info = Array{BiallelicVariantInfo,1}()
	open(bimfile) do bimf
		for line in each_trimmed_nonempty_line(bimf)
			chr,snpname,cm,pos,a1,a2 = getfields(line)
			push!(markers_info,BiallelicVariantInfo(chr,snpname,parse(Int32,cm),parse(Int32,pos),a1,a2,0,0))
		end
	end
	return markers_info
end

# With a dict we temporary relax constraints of equal number of colums per row
fam_dict() = Dict(
	:FID =>  data(Array{ASCIIString,1}()),
	:IID =>  data(Array{ASCIIString,1}()),
	:PaternalID =>  data(Array{ASCIIString,1}()),
	:MaternalID =>  data(Array{ASCIIString,1}()),
	:Sex => data(Array{Int8,1}()),
	:Phenotype => data(Array{Int8,1}()))

function fam_push!(fam::Dict, fields)
	push!(fam[:FID], shift!(fields))
	push!(fam[:IID], shift!(fields))
	push!(fam[:PaternalID], shift!(fields))
	push!(fam[:MaternalID], shift!(fields))
	push!(fam[:Sex], parse(Int8,shift!(fields)))
	push!(fam[:Phenotype], parse(Int8,shift!(fields)))
end

function fam_dataframe(dict::Dict) 
	df = DataFrame(dict) 
	for id in keys(dict)
		df[id] = map(i -> (i==-9 ? NA : i), df[id])
	end
	return df
end

function fam_readfile(famfile)
	fam=fam_dict() 
	open(famfile) do f
		for line in each_trimmed_nonempty_line(f)
			fields = getfields(line)
			@assert 6==length(fields) "$famfile malformed on line $line $famfile"
			fam_push!(fam,fields)
		end
	end
	return fam_dataframe(fam)
end

function readplink(p::PlinkFile) 
	markers_info = map_readfile(p.map)

	markers = Array{DataArray{BiallelicVariant,1},1}()
	fam=fam_dict()

	open(p.ped) do ped
		for line in each_trimmed_nonempty_line(ped) 
			fields = getfields(line)

			fam_push!(fam,fields)

			for i in 1:length(markers_info)
				a1=shift!(fields)
				a2=shift!(fields)

				if (a1!="0") 
					if markers_info[i].a1 == ""
						markers_info[i].a1 = a1	
					end

					if markers_info[i].a2 == ""  && (a1 != a2 || a1!=markers_info[i].a1)
						if markers_info[i].a1 == a1
							markers_info[i].a2 = a2
						else
							markers_info[i].a2 = a1
						end
					end
				end

				genoint = a1=="0" ? -9 : a1!=a2 ? 1 : a1==markers_info[i].a1 ? 0 : 2
				genotype = BiallelicVariant(markers_info[i],genoint)

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

	df = fam_dataframe(fam)

	for i in 1:length(markers)
		df[Symbol(markers_info[i].name)] = markers[i] 
	end

	return df
end

function byte_to_geno(byte) 
	map(0:3) do pos
		(byte >> (2*pos)) & 3
	end
end


function readplink(p::PlinkBinaryFile)
	markers_info = bim_readfile(p.bim)
	fam = fam_readfile(p.fam)
	markers = Array{DataArray{BiallelicVariant,1},1}()
	
	println(markers_info)

	# The first three bytes should be 0x6c, 0x1b, and 0x01 in that order. 
	bed=open(p.bed,"r")
	magic = readbytes(bed,3)

	@assert magic == [0x6c, 0x1b, 0x01] "magic header in .bed file does not look like a binary plink file"

	# The rest of the file is a sequence of V blocks of N/4 (rounded up) bytes each
	# where V is the number of variants and N is the number of samples. 
	# The first block corresponds to the first marker in the .bim file, etc.
	for mkr in markers_info
		last_byte_samples = nrow(fam) % 4 
		ablock = readbytes(bed,Integer(ceil(nrow(fam)/4)))
		four_sample_blocks = map(byte_to_geno,ablock)
		if last_byte_samples != 0
			for i in 1:(4-last_byte_samples)
				pop!(last(four_sample_blocks))
			end
		end

		genotypes=reduce(vcat,[],four_sample_blocks)
		
		push!(markers,DataArray{BiallelicVariant,1}[])
		for geno in genotypes
			if (geno == 0)
				push!(markers[end],BiallelicVariant(mkr,0))
			elseif (geno == 1)
				push!(markers[end],NA)
			elseif (geno == 2)
				push!(markers[end],BiallelicVariant(mkr,1))
			elseif (geno == 3)
				push!(markers[end],BiallelicVariant(mkr,2))
			end
		end
	end
	close(bed)

	for i in 1:length(markers)
		fam[Symbol(markers_info[i].name)] = markers[i] 
	end

	return fam
end

#dataframe(p::PlinkBinaryFile) = error("Not implemented yet")

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
