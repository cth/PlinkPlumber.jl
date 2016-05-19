#using PlinkPlumber
using Base.Test

# Make sure that plink is installed 
@test ismatch(r"PLINK.*",readall(`plink --version`))

function create_ped()
# Example from https://www.cog-genomics.org/plink2/formats
	test_ped = """
  	  1 1 0 0 1  0  G G  2 2  C C
  	  1 2 0 0 1  0  A A  0 0  A C
  	  1 3 1 2 1  2  0 0  1 2  A C
  	  2 1 0 0 1  0  A A  2 2  0 0
  	  2 2 0 0 1  2  A A  2 2  0 0
  	  2 3 1 2 1  2  A A  2 2  A A
	"""
 	 
	test_map = """
 	 1 snp1 0 1
 	 1 snp2 0 2
 	 1 snp3 0 3
	"""

	plinkstem=tempname()

	open(string(plinkstem,".ped"), "w") do f
		write(f,test_ped)
	end

	open(string(plinkstem,".map"), "w") do f
		write(f,test_map)
	end

	return plinkstem 
end

function create_bed(pedstem)
	plinkstem=tempname()
	run(`plink --file $pedstem --make-bed --out $plinkstem`) 
	plinkstem
end


function test_readplink()
	a_ped_file = create_ped()
	a_bed_file = create_bed(a_ped_file)

	pedread = readplink(plink(a_ped_file))
	bedread = readplink(plink(a_bed_file))

	println(pedread)
	println(bedread)
end

