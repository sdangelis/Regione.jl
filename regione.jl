### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ d582a3c1-c2e6-4519-a65e-950457169e31
using BED 

# ╔═╡ 7a478379-64f0-47b4-8724-33729f0a817e
using Distributions 

# ╔═╡ 257d3269-4586-4c97-b179-59565cd0c1f7
using GenomicFeatures

# ╔═╡ b158fecb-95bf-4a14-b06a-a5c74449c937
using Gtk

# ╔═╡ 9c1ad09e-5834-4528-8441-5e48106d3a8b
using Plots

# ╔═╡ 8af00452-364d-11ec-1f05-498c97e027e8
md"""
# Regione.jl

This notebook explores the possibility of a Julia drop-in replacement for the RegioneR Bioconductor package's overlap permuation test.

in particular, it aims to the implement a fast alternative to randomizeRegions
using the higher efficiency of Julia loop compared to R and the lower overheads 
introduced by GenomicRanges.jl over Bioconductor's GRanges


it might eventually implement the other 2 randomising functions

## Introduction
"""

# ╔═╡ fca01096-de01-45c6-8a0b-38b62ca22295
md"""
Load Dependencies
"""

# ╔═╡ 221bc89f-b5d3-428c-9fce-0e00d72cf679
tests = repeat([true], 16)

# ╔═╡ a5357f24-f3c5-4c38-b007-d30688f00839
md"""
## load beds 

We can use BED.jl to read bed files
This is the easy part
"""

# ╔═╡ 96af8850-364f-49e8-9f41-cad43d775ec9
features =  open(BED.Reader, open_dialog("Choose Feature File")) do reader
    IntervalCollection(reader, true)
end

# ╔═╡ bc98e879-1c37-4cf4-b882-e707f1485353
cnas = open(BED.Reader, open_dialog("Choose CNAs bed file")) do reader
    IntervalCollection(reader, true)
end

# ╔═╡ df7c40cd-4cd4-4aab-b174-13c0183fc51e
regions = open(BED.Reader, open_dialog("Chppse regions BED file")) do reader
    IntervalCollection(reader, true)
end

# ╔═╡ 4825b5a4-0acd-4e3a-8bf7-8caa143f5a65
md"""
NOTE: Generally sepeaking Loops in Julia are fast. Julia also supports vectorisation 

Vectorisation is not universally better and many base functions actually do the loops in Julia itself rather than vectorise and send the loop to C 
"""

# ╔═╡ c0203795-8ac9-420d-8271-19705ce23e53
md"""
## Functions to validate input 

Before we can continue, we need to find a way to check if two interval collections overlap with each other 

While I am sure there is a native implemetnation, for now we can use a naive appoach, 
simply builing on `Genomic.Features.is_overlapping` function

Indeed, there should be a way to intersect trees that improves on the O(a * b) cost of this implementation, with a and b being the number of intervals.
"""

# ╔═╡ 7112dc40-8fb8-47b9-84f6-df175a5f2f59
"""
extends GenomicFeatures.isoverlapping to  linearly searche if an interval A overlaps a collection b 
things are sorted so bisection search could be used if needed
things are also in a tree so maybe there's an even better way
"""
function GenomicFeatures.isoverlapping(a::GenomicFeatures.Interval{T}, b) where {T}
	for interval in b
		isoverlapping(a, interval) && return(true::Bool)
	end
	return(false::Bool)
end

# ╔═╡ b297f472-55b7-4a01-aa60-da1e993588c3
"""
Linearly check if all intervals in collection a overlap with collection b
Here the double linearity becomes O(a*b)
Hower interval collections seem to have little if any overheads compar
"""
function alloverlapping(a, b)
	for interval in a
		!isoverlapping(interval, b) && return(false::Bool)
	end
	return(true::Bool)
end

# ╔═╡ 9bead016-407f-4d5f-8912-7a5132241ec6
md"""
Let's test this on synthetic intervals
"""

# ╔═╡ fabdc60e-0365-4cf6-ac0b-cdd73d56a0ea
begin
test1 = GenomicFeatures.IntervalCollection([GenomicFeatures.Interval("chr1", 1, 100),GenomicFeatures.Interval("chr1", 200, 300),GenomicFeatures.Interval("chr1", 400, 500)]) # interval b
test2 = GenomicFeatures.Interval("chr1", 1150, 1180) # False
test3 = GenomicFeatures.Interval("chr1", 250,280) # true
test4 =  GenomicFeatures.IntervalCollection([GenomicFeatures.Interval("chr1", 50, 80),GenomicFeatures.Interval("chr1", 250, 280)]) # true
test5 =  GenomicFeatures.IntervalCollection([GenomicFeatures.Interval("chr1", 50, 80),GenomicFeatures.Interval("chr1", 2250, 2280)]) # false (some true) 
test6 = GenomicFeatures.IntervalCollection([GenomicFeatures.Interval("chr1", 2150, 2180),GenomicFeatures.Interval("chr1", 2250, 2280)]) # false
test7 =  GenomicFeatures.Interval("chr1", 99,280) # false
end

# ╔═╡ 07e85de4-3e17-45be-a8ea-8b3d200cca7a
tests[1] = ([false, true, true] == map(x->isoverlapping(x,test1), [test2, test3, test7]))

# ╔═╡ 47e62456-1990-4f99-a23e-7f34781e6450
tests[2] == ([true, false, false]  == map(x->alloverlapping(x,test1), [test4, test5, test6]))

# ╔═╡ d99cfcf1-202a-4b1a-893a-0e4eb6eddb94
tests[3] = (true == alloverlapping(cnas, cnas))
# any collection should always overlap all with itself

# ╔═╡ 5f11206b-9689-4998-9dec-8a967db53bd6
md"""
At this point of the story we have a functioning way to see if all intervals in a collection overlap with any interval in another collection, **IGNORING STRAND but NOT SEQUENCES**. 
However, our beds come in a collection with mutliple sequences.

We can subset for the sequence we want, but while iterable, this type is technically distinct from an interval collection and we will need to
edit the type information on the function. In theory we should be able to make the functions generic and let Julia deal with it.

But for simplicity we can explicitly transform an interval tree into interval collection by iteration or in one go which should be faster
Fun Fact! does not seem to be the case

Again Julia's loop optimisations means loops are no longer the enemy 

Interestingly this operation seems to be fast enough to run without worrying about changing all functions to deal with this type
"""

# ╔═╡ 91745a67-8b05-48b5-b42a-cc5893fef18b
function vec_getcollection(collection, sequence) 			GenomicFeatures.IntervalCollection([GenomicFeatures.Interval(i.seqname, 
				i.first, i.last, i.strand, i.metadata ) for i in collection.trees[sequence]])
end

# ╔═╡ fd60dc2e-af8f-48ba-a2da-ab8e62ad3c83
function iter_getcollection(collection::IntervalCollection{T}, sequence) where {T}
tmp = IntervalCollection{T}()
	for i in collection.trees[sequence]
		push!(tmp, i)
	end
return(tmp)
end

# ╔═╡ 7e78757e-464a-4d22-9adc-45fa1b799ff2
test11 = vec_getcollection(cnas, "chr1")

# ╔═╡ 34252bfe-3b85-4795-9ed0-1c68a46bc2f2
test12 =  iter_getcollection(cnas, "chr1")

# ╔═╡ 1657a4c7-e489-45a4-b52c-b7bb12bb7c63
tests[4] == (test11 == test12)

# ╔═╡ b10ce41f-d915-405a-b681-321ffdc48c12
md"""
everthing looks good, oddly both functions are comparable, with the loop being very fast. Likely thanks to in place operations. From a R/Python perspective loops are no longer the enemy

## Overlap Counting

Now we can subset ranges and calculate if all overlap, it is time to think about counting overlaps

Here if we follow the naive implementation we'd get a O(a*b) again. 
"""

# ╔═╡ a28e9dce-e2ce-4f07-a154-10cf5ee40326
"""
Linearly count how many intervals in collection a overlap with collection b.\n
Note countoverlapping(a,b) != countoverlapping(b,a) as one overlap in a can actually overlap with >1 intervals in b
"""
function countoverlapping(a::GenomicFeatures.IntervalCollection{T} ,b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	overlaps::Int64 = 0
	for interval in a
		isoverlapping(interval, b) && (overlaps += 1) 
	end
	return(overlaps)
end

# ╔═╡ d3510b47-8b89-4e3e-b8f5-de5b68addd72
md"""
we can borrow the test set 1 to check this function does what we expect
"""

# ╔═╡ eccaee6d-3534-4daa-aaff-3d6023a5b189
tests[5] = (true == (cnas.length == countoverlapping(cnas,cnas)))
# overlap of itself should return length 

# ╔═╡ 4fd1cf67-6fc5-4045-9a37-23dd3c73ed01
md"""
## In 

However, ovelapping by itself does not tell us whether a collection a is a subset of collection b as overlapping only guarantees that interval a either start or ends within interval b, not both. We can check this naively and linearly with a series of isin functions. 

Julia's multiple dispatch system is kind enough to let us safely override the built-in `in` operator with our own, enabling us to do quite powerful things

Here the fact that intervals are sorted really could have given us more effective methiods to do this rather than again linear search 

Nonetheless, we might get away with this.
"""

# ╔═╡ df44379b-40f1-4046-8039-a653c8b79e7e
"""
checks if interval a is fully contained in interval b
"""
function Base.in(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S,T}
	return (a.first >= b.first && a.last <= b.last)
end

# ╔═╡ 00a5c4d3-768c-4904-85ec-da9eb42cb1ea
"""
linearly checks if interval a is contained in intervalcollection b 
"""
function Base.in(a::GenomicFeatures.Interval{S},b::GenomicFeatures.IntervalCollection{T}) where {S,T}
	for interval in b
		in(a, interval) && return(true::Bool)
	end
	return(false::Bool)
end

# ╔═╡ f93f9961-4d63-4c27-b6e1-c1127da7498c
"""
linearly checks if all intervals in a are fully contained in intervalcollection b 
"""
function Base.in(a::GenomicFeatures.IntervalCollection{S},b::GenomicFeatures.IntervalCollection{T}) where {S,T}
	for interval in a
		!in(interval, b) && return(false::Bool)
	end
	return(true::Bool)
end

# ╔═╡ b5e5d87a-ad5c-491e-94b7-1beb865b925d
# check if one interval is in the other 
tests[6] =  ([true,false] == [in(test3,test7),in(test2,test3)])

# ╔═╡ 528f5030-717e-4307-9be3-4bce601563d8
# check if one interval is in the collection 
tests[7] == ([false, false, true] == [in(test7,test1), in(test2,test1), in(test3,test1)])

# ╔═╡ 1b74fef9-c507-4b01-9143-6a944dd566af
# check if one collection is in the other
tests[8] = ([true, false, false] == [in(test4,test1), in(test5, test1), in(test6, test1)])

# ╔═╡ e7885943-2e64-4b32-af85-0394feb0290e
tests[9] = (false == (test2 in test1))

# ╔═╡ c8fa0127-34bf-4900-b40c-446a464b50a1
md"""
All good here, 

## Generating Random Regions

Now we need to start generating the random regions.
To do so we can harness Distributions.jl ability to generate statitscal distributions on the fly.
Given we want to support custom bed regions for the randomisation, these will necessarily contain gaps. If we allow features to go in these gaps
Briefly we can model the start position as a Mixture model with a categorical distribution of random uniform distributions. The categorical distribution will allow
us to sample from either distribution.

For more information on where this idea comes from see [this Stack Overflow thread](https://stackoverflow.com/questions/49924115/julia-random-number-in-different-intervals/49924852#49924852)

Here comes another interesting question. What are we randomising over?
Are we simply interseted in the overlap within
Or are we interested in the postion of CNAs across chromosomes

Likely the answer will vary by the circumstance.

Another question is what heruistic to use to make sure that we move intervals of variable length in the valid regions. i.e we don't get an interval of the right strat but that ends past the valid regions

there are 2 contrasting approaches here 

1. sample start location until the end 
2. generate an array of distributions for each given length and store in a key-value pair, using each distribution as needed.

Which approach to take depends on the relative speed of generating/accessing multiple distributions vs checking and resampling. Intuitively, for larger features the second should be better while the first should work for smaller features.
"""

# ╔═╡ f8f1714f-0ea2-4779-a1e2-2db68f992152
md"""
## Generate Distribution

generating a discrete mixture model is the easy step. Straight from Stack Overflow 

To check we are doing what we think we are we can create a mixture model from 1 to 10 based on two uniform distrbutions [1-4, 5-10] with relative priors [0.4,0.6]

Given there is no gaps and the categorial distributions match length we should get something equivalent to a uniform 1 to 10
"""

# ╔═╡ 89769f66-9cf8-4619-ab57-256bd683f695
begin 
test31 =  MixtureModel([DiscreteUniform(1,4), DiscreteUniform(5,10)], Categorical([0.4,0.6]))
test32 = DiscreteUniform(1,10)
end

# ╔═╡ 9077f8e3-7dcc-46dc-85c8-6ffd7c9c5c0b
test33 = (mean(test31) == mean(test32)) == (std(test31) == std(test32))

# ╔═╡ 4a249ad1-cb44-4f9c-b61c-f560e350824a
test34 = isapprox(pdf(test31), pdf(test32)) # Note we need isapprox due to floats

# ╔═╡ 97d9b1b4-453b-4e2b-a7dd-4dc941900c00
tests[10] = (test33 == test34 == true) 

# ╔═╡ c13c7337-db5d-4177-8a6c-79ff3c1173f2
pdf(test31) == pdf(test32)

# ╔═╡ db80c34b-7970-4d72-89b4-c90091c1f553
pdf(test31), pdf(test32) # not the same, yet look the same. Mystery of the Floats 🤷 

# ╔═╡ 03bae442-29ee-4272-b630-d051956a60ab
isapprox(pdf(test31), pdf(test32)) #A-Ha!, Hello Floats my old friends

# ╔═╡ 869ec51c-e6f5-4733-84f8-6b006dc0f294
md"""
### Approach 1

same distribution with resampling until we find an interval that fits. Here we rely on our custom in method.

To allow the whole system to fail safely, should we get to the point where we can generate random regions that cannot belong to the regions collection we need to allow only a finte number of tries before giving an error

But how often we expect this to happen? For instnace if the length of an interval to randomise exceeds the length of the largest interval in the allowed regions, or if we supply intervals in a different sequence that does not exist in the regions, or we have a mismatched sequence/distribution 
"""

# ╔═╡ 8d99a321-f9c6-4da9-921d-55bf50034dc5
"""
Generates a mixture model from a given genomic inteval collection
"""
function generatedistribution(collection)
d = Vector{Distributions.DiscreteUniform}(undef, 0)
l = Vector{Int64}(undef,0)
for i in collection
	push!(d, DiscreteUniform(i.first, i.last))
	push!(l, i.last-i.first)
end
l = l/sum(l)
return MixtureModel(d, Categorical(l))
end

# ╔═╡ 01aae5a0-f208-477d-974b-3945f0d56fd0
M = generatedistribution(test1)

# ╔═╡ ab81fd7e-f6fe-4c25-b5b1-9f19f4d9cddb
test6

# ╔═╡ 245afe84-62df-4496-b8ff-c4d586fb0f47
"""
generates a new random interval from a given distribution of possible start locations ensuring it is contained in the given regions

Optionally intervals in the collection can be disallowed with allow_overlap = false, default: true. 
"""
function _randominterval(interval, distribution, regions; collection=nothing, allow_overlap = true) 
il = interval.last - interval.first
new = GenomicFeatures.Interval
max_iter = 1000
for i in 1:max_iter
	newfirst = rand(distribution,1)
	new = GenomicFeatures.Interval(interval.seqname, newfirst[1], newfirst[1]+il,
			interval.strand, interval.metadata)
	if allow_overlap
			new in regions && return new
	else 
			(!(isoverlapping(new, collection)) && new in regions) && return new
	end
end 
error("Error: can't randomise regions with $max_iter tries,
		Please consider increasing max_iter or use randomisation strategy 2,
		Please also check the intervals come from the same sequnce used to generate the distribution from the regions")
end

# ╔═╡ db717a7f-bc68-47fb-b894-ff97bd1bd1ec
begin 
	test41 = _randominterval(test2,M, test1) 
	tests[11] = ((true, false) == (test41 in test1, test41 == test2))
end

# ╔═╡ 92248a65-e3fe-4071-b410-ea853b4c4c11
md"""
so far we have a way to generate a random interval. we can now extend this (with a loop) to generate a new collection altogether
"""

# ╔═╡ bc322aa5-c57f-487c-9e5a-c907e3582ef7
"""
Generates randomised regions from an interval collection according to a start postion  distribution. It guarantees that the new intervals are contained in the given regions.

Note for performance reasons this functions assumes that the interval collection contains intervals belonging to a single sequence that is the same sequence used for the distribution. The distribution should be derived from the regions, else the randomiser will (safely) return an error 

New intervals overlapping with other intervals in the nascient collection are disallowed by default. To allow overlaps use `allow_overlap=true`
"""
function randomiseregions(collection::IntervalCollection{T}, distribution, regions; allow_overlap = false) where {T}
	new_collection = GenomicFeatures.IntervalCollection{T}()
	for i in collection
		push!(new_collection, _randominterval(i, distribution, regions, collection = new_collection, allow_overlap = allow_overlap))
	end 
	return new_collection
end 

# ╔═╡ 1621bce6-e8bb-4d70-8784-21aa80474af3
test1

# ╔═╡ 1dc3eec5-fe44-474a-a8df-007102a522d3
randomiseregions(test4, M, test1)

# ╔═╡ dadb9d11-80ba-4e85-b49d-e7671672dc6e
(randomiseregions(test4, M, test6))
# looks like this is now safe enough for use!

# ╔═╡ 983d0f16-1484-44c1-818d-5827ab582e5f
md"""
### Randomise at scale

Let's now test this all with the larger intervals from the bed file.
"""

# ╔═╡ 5abe2ecd-51bf-4854-b287-1ce85ab30347
features_chr1 = iter_getcollection(features, "chr1")

# ╔═╡ 192c4600-a9c1-4dbb-a74c-6352bb125059
regions_chr1 =  iter_getcollection(regions, "chr1")

# ╔═╡ 12afb103-311c-4e2b-8243-1ca18f7dd14a
tests[16] = (countoverlapping(features_chr1, regions_chr1) != 
countoverlapping(regions_chr1, features_chr1))
# This should really NOR be true

# ╔═╡ 81be70da-010d-45b6-9e8e-d5500dcbb867
features_chr1 in regions #check if we have all features

# ╔═╡ c80a9688-3ceb-45e2-8c18-93d0599f03fc
features in regions_chr1 # as a test, this should return false

# ╔═╡ bbdecd68-45cc-4894-8d19-cd13a51bdedc
cnas_chr1 = iter_getcollection(regions, "chr1")

# ╔═╡ bd33b46b-bfdc-449e-a614-cc04688f5b89
begin 
	o = countoverlapping(features_chr1, cnas)
	D = generatedistribution(regions_chr1)
	o2 = Vector{Int64}(undef,1)
	for i in 1:1
		random = randomiseregions(features_chr1, D, regions_chr1)
		o2[i] = countoverlapping(random, cnas)
	end
o2
end

# ╔═╡ 58e67962-6cec-4d01-b4f8-ce003c0e4967
tmp_cnas = iter_getcollection(cnas, "chr1")

# ╔═╡ 9cdcd137-a2fb-4f82-98a3-e613e1a80f54
random = randomiseregions(tmp_cnas,	D, regions_chr1) # now safe enough to throw an error! but why is this not working? likely length of the segment is too big for the distribution. Gotta check this in the final function

# ╔═╡ c96c018b-d4a4-4ebf-925c-bbf9bcb4fc8f
tmp_cnas in regions_chr1

# ╔═╡ f6b3bb75-3565-41f0-b592-81b2b36eb651
cnas in regions

# ╔═╡ 35f82e76-05c0-4aa4-823e-16187d90ab25
md"""
This is the performance of the randomiser per 1 sequence.
The sharper eyed will have noted we have not checked for overlaps and are artificially limiting this to one chromosome. It should be easier to define all our functions this way and then let the permutation test function deal with it.

Morevoer, at the moment this cannot randomise long CNAs effectively as there are many random intervals that need to be rejected due to lenght. This can slow things down a lot.
"""

# ╔═╡ 9de81723-c470-4e53-8fb8-ea2c1a889c1f
md"""
## Permutation Test 

The next step now is to package what we've done so far into a permutation test.

For now we can borrow RegioneR's P value calculations but in the future I would like to know where this comes from. It seems too simple?

This value is derived as ``p = \frac{n+1}{i+1}``, where n is the number of random samples that are higher or lower (based on alternative) than the observed value and i is the total number of observations.
"""

# ╔═╡ 2b1689f4-89a7-42a3-916d-5731dc91b5b3
ntimes::Int64= 1e4 # Define number of iterations

# ╔═╡ d5f7057c-da96-405e-9ad8-a8032c7c47a1
Text("The Minimum P-Value with $ntimes iteration is $(1/(ntimes+1))") # get max p-value

# ╔═╡ 926897bd-edca-4961-975b-7117b7d000e9
Dict([(1, 1),(2, 2)])

# ╔═╡ 07569676-36be-4b54-8105-c28e97bb872a
"""
Calculates the p value for the alternative hypothesis that the observed value is either more or less than what we expect from the Vecor of random values.
If no alternative is provided, an alternative will be automatically selected 
on the basis of the data. 
Please see regioneR for full details
"""
function calculatep(obs, rand; alternative = "Auto")
	if alternative == "Less"
		p = (length(filter(x-> x <= obs, rand)) + 
          1)/(length(rand) + 1)
	elseif alternative == "More" 
		p = ((length(filter(x-> x >= obs, rand)) + 
          1)/(length(rand) + 1))
	else
		obs <= mean(rand) && return(calculatep(obs, rand; alternative = "Less"))
		return(calculatep(obs, rand, alternative ="More"))
	end
	return(p = p, alternative = alternative)
end 

# ╔═╡ 122767a3-3348-4b4c-990f-315dfe7c30ca
calculatep(2,[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,11,1,1,1,1,1,1,1]; alternative = "Auto")

# ╔═╡ 156e04fd-cf99-4e11-97a9-eb334a85a2cb
begin 
test51 = calculatep(2,repeat([1],9))
test52 = calculatep(2,repeat([3],9))
test53 = calculatep(2,repeat([2],9))
test54 = calculatep(2,repeat([1],1000))
end

# ╔═╡ 636135ce-bfd0-4786-ad99-25af1832fe2f
test51,test52,test53,test54

# ╔═╡ ea7dd844-8d79-4d47-98ce-c91fd3fa04b6
begin 
tests[12] =  test51 == (p = 0.1, alternative = "More")
tests[13] = test52 == (p = 0.1, alternative = "Less")
tests[14] = test53 == (p = 1.0, alternative = "Less")
tests[15] = test54 == (p = 0.000999000999000999, alternative = "More")
end 


# ╔═╡ 55317d1e-88de-4b61-8b1a-616463ab7299
md"""
We now have a way to calculate the P value of our overlap permutation test and what is the minimum value it can take.

For now we are not considering the Z test as I really am not sure we can do that with 1 observation (and HypothesisTesting.jl doesn't seem to like that) 

Alternatively, the KS test seems to be promising, as we can use to to ask wether a vector comes from a given distribution. Again is this appropriate for 1 sample? 
"""

# ╔═╡ a4a7ecb8-293b-47a4-a45a-b0fbc7f9fc7e
md"""
## Putting Everything Together 

Now that we have all the elements ready we can put everything in our test together,

Note:
1. for now we are only supporting randomisation strategy 1. 
2. we are not implementing the opportunity to pass a custom evaluation or randomisation function.

Both should be fairly trivial to do 
"""

# ╔═╡ 0bf2dc1c-084f-4714-8bb6-11e631767199
"""
Simple macro to get the name of a variable, thanks to [Micheal Ohlrogge (https://stackoverflow.com/users/3541976/michael-ohlrogge)
"""
macro Name(arg)
	string(arg)
end

# ╔═╡ 213e99b1-13ba-48cb-9838-f80077d85992
"""
A simple structre and constructor to store results of the overlap permuation test
"""
struct OverlapTestResult
	iterations::Int64
	obs::Int64
	rand::Vector{Int64}
	p_val::Float64
	alternative::String
	min_p_val::Float64
	tested_regions::Union{String, Nothing}
	randomised_regions::Union{String, Nothing}
	function OverlapTestResult(iterations, obs, rand, p_val, alternative, names)
		new(iterations, obs, rand, p_val, alternative, (1/(iterations+1)), names[1], names[2])
	end
end

# ╔═╡ d2773d45-b995-48b4-bcf8-754341dc3dcd
"""
TBD (to be documented 
"""
function overlaptest(col1, col2, regions, iterations; bed_names = (nothing,nothing), allow_overlap=false)
	name =
	# count base overlaps 
	obs = countoverlapping(col1,col2)
	if !(col2 in regions)  
		@warn "Intervals to randomise are not a subset of target regions for the randomisation"
	end
	# iterate per chromosome, generating all distributions 
	D = Dict{String, Distributions.MixtureModel}()
	for seq in keys(regions.trees)
		push!(D, seq => generatedistribution(regions.trees[seq]))
	end
	
	# iterate per chromosome get all regions 
	
	r = Dict{String, 
IntervalCollection{}}()
	for seq in keys(regions.trees)
		push!(r, seq => iter_getcollection(regions, seq))
	end
	
	a = Dict{String, 
IntervalCollection{}}()
	for seq in keys(col1.trees)
		push!(a, seq => iter_getcollection(col1, seq))
	end
	
	b = Dict{String,
IntervalCollection{}}()
	for seq in keys(col2.trees)
		push!(b, seq => iter_getcollection(col2, seq))
	end
	 

	# over the number of iterations
	rand = Array{Int64}(undef, iterations)

	for i in 1:iterations
		tmprand = Array{Int64}(undef,0)
		for seq in keys(col2.trees)
			# run chromosome by chromosome.
			rand_b = randomiseregions(b[seq],D[seq],r[seq]; allow_overlap = allow_overlap)
			push!(tmprand, countoverlapping(a[seq], rand_b))
		end 
		rand[i] = sum(tmprand)
	end 
	# get p value for overlap test
	p_val = calculatep(obs, rand)
		
	return OverlapTestResult(iterations, obs, rand, p_val[:p], p_val[:alternative], bed_names)
end

# ╔═╡ ad624619-cc2e-42c5-ae62-817771f600ee
md"""
Step 1: check if we pass unit tests.
let's collect the results and compare
"""

# ╔═╡ 951cde36-9000-4d8f-9cb5-81401be51671
pass = unique(tests)

# ╔═╡ 94daad38-215e-4f54-8165-af5f8a373804
md"""
If all is good we can finally run our overlap permuatation test function on synthetic data 
"""

# ╔═╡ aa3c015d-371a-444d-9b9c-98d2aef53d1e
md"""
## All Systems go 🚀

Pluto is a really amazing interactive notebook, so by now the state of this function depends on everything done so far and will get dynamically updated every time we change anything above.

Note this will run at every modification, even if something is broken. Fortunately we did unit tests as we went along so if something breaks, wan use the results of our tests to check all is good before running the final overalp premutation step
"""

# ╔═╡ 990a8734-a7a9-4d11-b10f-447912e47e60
it_go = all(tests .== true)

# ╔═╡ fe573923-fd00-4c76-ac58-eedac559b6a6
md"""
### Integration test

Before we do the test on real data, we can check everything is what we expect
"""

# ╔═╡ 38038448-d16d-4e0c-8ac1-233aca2635e3
test61 = IntervalCollection([GenomicFeatures.Interval("chr1", 1, 500000)])

# ╔═╡ c0a4b035-baf1-452f-a608-95c4eeaa6549
it_go && (test62 = overlaptest(test4, test5, test61, 100000))

# ╔═╡ d7666993-dc60-447c-ab90-696b2c9862c1
test63 = test62.p_val < 0.005 # test 62 should  reliably return a p-value below 0.005

# ╔═╡ 649d0694-704f-4de5-a882-338a276968c8
test64 = overlaptest(test1,test1, test61, 10000;allow_overlap = false) 

# ╔═╡ 9fabd4e0-b733-48e2-81bf-b0c6b5a4dae1
test65 = test64.p_val == test64.min_p_val 	# overlap with itself should be as significant as it gets 

# ╔═╡ 1945ba42-eb1d-469f-bf53-c5b8e7573b54
go =  all([test63,test65] .== true)

# ╔═╡ 4a4d0612-779f-4715-beb6-913f8023cccd
md"""
## Showtime

with all unit and integration tests working, we can now apply the test to real data.

The bed file CNAs here contains all regions called, labelled as CNV 1000. Hence it will overlap all large genes. However, as the gene position will get reshuffled the number of segments 
that overlap at least 1 gene will change as genes get moved  
"""

# ╔═╡ 5f46bdeb-31af-417e-9e80-291898252734
go && (result1 = overlaptest(cnas, features, regions, 10000; allow_overlap = false))

# ╔═╡ 7c2560f8-c8cc-4658-ab17-77d8322b2dfd
begin 
histogram(result1.rand, bins = :scott, color_palette=["grey75"]) # Ok, not the prettiest but it can fly for now.
end

# ╔═╡ a2b2a498-524c-47f2-a15f-bdc0768f820a
md"""
### Performance 

at 100,000 iterations per 2 minutes($(100000/108) iteration per second) we compare quite favourably with regioneR despite not including plotting or Z test 

In a similar context regioneR does ($(100/180) iteration per second), plotting + Z test included (but they should be fairly lightweigth)


Hovever this has a cost, namely safety/correctness. While RegioneR builds upon GRanges and has an extensive way to validate intervals, sequences, automatically extract sequence names and esure they are comparible with the sepcifie genome, this notebook assumes most are dealt with manually.

it realistically only throws an error where the situation would lead to an infinite loop. For instance, there is nothing stopping you from using hg38 regions to randomise hg19! Although the opposite would possibly trigger an error as there is no way to randomise something on a chromosome that does not exist in the regions

Moreover, I have likely cut multiple corners in a way that will affect the results here so this is not an autoritative benchmark

**Use at your own risk!**

### what about randomising regions 

Can we do the opposite and actually randomise the regions, with the expectation that the number of overlaps is constant?

Well turns out the bed file contains some segments spanning the centromere so at the moment we expect the test to error. 
"""

# ╔═╡ 95eaff97-3e79-4b3b-8265-0bd77bd31183
result2 = overlaptest(features, cnas, regions, 100; allow_overlap = false)
# this is at the moment impossible and rightly errors.

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BED = "8e4a8c10-cb6b-11e8-08d2-83478d609d67"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
GenomicFeatures = "899a7d2d-5c61-547b-bef9-6698a8d05446"
Gtk = "4c0ca9eb-093a-5379-98c5-f87ac0bbbf44"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
BED = "~0.2.1"
Distributions = "~0.25.22"
GenomicFeatures = "~2.0.4"
Gtk = "~1.1.9"
Plots = "~1.23.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ATK_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "58c36d8a1beeb12d63921bcfaa674baf30a1140e"
uuid = "7b86fcea-f67b-53e1-809c-8f1719c154e8"
version = "2.36.1+0"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[BED]]
deps = ["Automa", "BGZFStreams", "BioGenerics", "ColorTypes", "FixedPointNumbers", "GenomicFeatures", "Indexes", "TranscodingStreams"]
git-tree-sha1 = "b353d994d389c3b2dc64ce83a5e8f5213f1e19e5"
uuid = "8e4a8c10-cb6b-11e8-08d2-83478d609d67"
version = "0.2.1"

[[BGZFStreams]]
deps = ["CodecZlib", "Test"]
git-tree-sha1 = "b0e322aef9f1895a09a0c48a118111fd509f666a"
uuid = "28d598bf-9b8f-59f1-b38c-5a06b4a0f5e6"
version = "0.3.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BioGenerics]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "6d3f3b474b3df2e83dc67ad12ec63aee4eb5241b"
uuid = "47718e42-2ac5-11e9-14af-e5595289c2ea"
version = "0.1.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "97f1325c10bd02b1cc1882e9c2bf6407ba630ace"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.12.16+3"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3fcfb6b34ea303642aee8f85234a0dcd0dc5ce73"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.22"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cafe0823979a5c9bff86224b3b8de29ea5a44b2e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.61.0+0"

[[GTK3_jll]]
deps = ["ATK_jll", "Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Libepoxy_jll", "Pango_jll", "Pkg", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXcomposite_jll", "Xorg_libXcursor_jll", "Xorg_libXdamage_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "Xorg_libXrender_jll", "at_spi2_atk_jll", "gdk_pixbuf_jll", "iso_codes_jll", "xkbcommon_jll"]
git-tree-sha1 = "f2922cb0105ba076b560e3922910f8e3749aa91a"
uuid = "77ec8976-b24b-556a-a1bf-49a033a670a6"
version = "3.24.30+0"

[[GenomicFeatures]]
deps = ["BioGenerics", "DataStructures", "IntervalTrees"]
git-tree-sha1 = "9f3f122571d4b279e5176eded5ece804fde93694"
uuid = "899a7d2d-5c61-547b-bef9-6698a8d05446"
version = "2.0.4"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[Gtk]]
deps = ["Cairo", "Cairo_jll", "Dates", "GTK3_jll", "Glib_jll", "Graphics", "Libdl", "Pkg", "Reexport", "Serialization", "Test", "Xorg_xkeyboard_config_jll", "adwaita_icon_theme_jll", "gdk_pixbuf_jll", "hicolor_icon_theme_jll"]
git-tree-sha1 = "d679fc90e75984e7b8eb25702eff7cf4718a9c10"
uuid = "4c0ca9eb-093a-5379-98c5-f87ac0bbbf44"
version = "1.1.9"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[Indexes]]
deps = ["BGZFStreams", "BioGenerics", "GenomicFeatures", "TranscodingStreams"]
git-tree-sha1 = "275bce824b40fd2e70358e0a652ba1b34172f240"
uuid = "4ffb77ac-cb80-11e8-1b35-4b78cc642f6d"
version = "0.1.3"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalTrees]]
deps = ["InteractiveUtils", "Profile", "Random", "Test"]
git-tree-sha1 = "6c9fcd87677231ae293f6806fad928c216ab6658"
uuid = "524e6230-43b7-53ae-be76-1e9e4d08d11b"
version = "1.0.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libepoxy_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "ae76616c0e0cbcbedc68d0121fe409598efc0b05"
uuid = "42c93a91-0102-5b3f-8f9d-e41de60ac950"
version = "1.5.8+0"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9bc1871464b12ed19297fbc56c4fb4ba84988b0d"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.47.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "f19e978f81eca5fd7620650d7dbea58f825802ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "ca7d534a27b1c279f05cd094196cb70c35e3d892"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.23.2"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMD]]
git-tree-sha1 = "9ba33637b24341aba594a2783a502760aa0bff04"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.3.1"

[[ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "9cc2955f2a254b18be655a4ee70bc4031b2b189e"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "95072ef1a22b057b1e80f73c2a89ad238ae4cfff"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.12"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcomposite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll"]
git-tree-sha1 = "7c688ca9c957837539bbe1c53629bb871025e423"
uuid = "3c9796d7-64a0-5134-86ad-79f8eb684845"
version = "0.4.5+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdamage_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll"]
git-tree-sha1 = "fe4ffb2024ba3eddc862c6e1d70e2b070cd1c2bf"
uuid = "0aeada51-83db-5f97-b67e-184615cfc6f6"
version = "1.1.5+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libXtst_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll"]
git-tree-sha1 = "0c0a60851f44add2a64069ddf213e941c30ed93c"
uuid = "b6f176f1-7aea-5357-ad67-1d3e565ea1c6"
version = "1.2.3+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[adwaita_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "hicolor_icon_theme_jll"]
git-tree-sha1 = "37c9a36ccb876e02876c8a654f1b2e8c1b443a78"
uuid = "b437f822-2cd6-5e08-a15c-8bac984d38ee"
version = "3.33.92+5"

[[at_spi2_atk_jll]]
deps = ["ATK_jll", "Artifacts", "JLLWrappers", "Libdl", "Pkg", "XML2_jll", "Xorg_libX11_jll", "at_spi2_core_jll"]
git-tree-sha1 = "f16ae690aca4761f33d2cb338ee9899e541f5eae"
uuid = "de012916-1e3f-58c2-8f29-df3ef51d412d"
version = "2.34.1+4"

[[at_spi2_core_jll]]
deps = ["Artifacts", "Dbus_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXtst_jll"]
git-tree-sha1 = "d2d540cd145f2b2933614649c029d222fe125188"
uuid = "0fc3237b-ac94-5853-b45c-d43d59a06200"
version = "2.34.0+4"

[[gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "0facfc4bfd873c21b83a053bbf182b9ef19c69d8"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.6+0"

[[hicolor_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b458a6f6fc2b1a8ca74ed63852e4eaf43fb9f5ea"
uuid = "059c91fe-1bad-52ad-bddd-f7b78713c282"
version = "0.17.0+3"

[[iso_codes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5ee24c3ae30e006117ec2da5ea50f2ce457c019a"
uuid = "bf975903-5238-5d20-8243-bc370bc1e7e5"
version = "4.3.0+4"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─8af00452-364d-11ec-1f05-498c97e027e8
# ╟─fca01096-de01-45c6-8a0b-38b62ca22295
# ╠═d582a3c1-c2e6-4519-a65e-950457169e31
# ╠═7a478379-64f0-47b4-8724-33729f0a817e
# ╠═257d3269-4586-4c97-b179-59565cd0c1f7
# ╠═b158fecb-95bf-4a14-b06a-a5c74449c937
# ╠═9c1ad09e-5834-4528-8441-5e48106d3a8b
# ╠═221bc89f-b5d3-428c-9fce-0e00d72cf679
# ╟─a5357f24-f3c5-4c38-b007-d30688f00839
# ╠═96af8850-364f-49e8-9f41-cad43d775ec9
# ╠═bc98e879-1c37-4cf4-b882-e707f1485353
# ╠═df7c40cd-4cd4-4aab-b174-13c0183fc51e
# ╟─4825b5a4-0acd-4e3a-8bf7-8caa143f5a65
# ╟─c0203795-8ac9-420d-8271-19705ce23e53
# ╠═7112dc40-8fb8-47b9-84f6-df175a5f2f59
# ╠═b297f472-55b7-4a01-aa60-da1e993588c3
# ╟─9bead016-407f-4d5f-8912-7a5132241ec6
# ╠═fabdc60e-0365-4cf6-ac0b-cdd73d56a0ea
# ╠═07e85de4-3e17-45be-a8ea-8b3d200cca7a
# ╠═47e62456-1990-4f99-a23e-7f34781e6450
# ╠═d99cfcf1-202a-4b1a-893a-0e4eb6eddb94
# ╟─5f11206b-9689-4998-9dec-8a967db53bd6
# ╠═91745a67-8b05-48b5-b42a-cc5893fef18b
# ╠═fd60dc2e-af8f-48ba-a2da-ab8e62ad3c83
# ╠═7e78757e-464a-4d22-9adc-45fa1b799ff2
# ╠═34252bfe-3b85-4795-9ed0-1c68a46bc2f2
# ╠═1657a4c7-e489-45a4-b52c-b7bb12bb7c63
# ╟─b10ce41f-d915-405a-b681-321ffdc48c12
# ╠═a28e9dce-e2ce-4f07-a154-10cf5ee40326
# ╟─d3510b47-8b89-4e3e-b8f5-de5b68addd72
# ╠═eccaee6d-3534-4daa-aaff-3d6023a5b189
# ╠═12afb103-311c-4e2b-8243-1ca18f7dd14a
# ╟─4fd1cf67-6fc5-4045-9a37-23dd3c73ed01
# ╠═df44379b-40f1-4046-8039-a653c8b79e7e
# ╠═00a5c4d3-768c-4904-85ec-da9eb42cb1ea
# ╠═f93f9961-4d63-4c27-b6e1-c1127da7498c
# ╠═b5e5d87a-ad5c-491e-94b7-1beb865b925d
# ╠═528f5030-717e-4307-9be3-4bce601563d8
# ╠═1b74fef9-c507-4b01-9143-6a944dd566af
# ╠═e7885943-2e64-4b32-af85-0394feb0290e
# ╟─c8fa0127-34bf-4900-b40c-446a464b50a1
# ╟─f8f1714f-0ea2-4779-a1e2-2db68f992152
# ╠═89769f66-9cf8-4619-ab57-256bd683f695
# ╠═9077f8e3-7dcc-46dc-85c8-6ffd7c9c5c0b
# ╠═4a249ad1-cb44-4f9c-b61c-f560e350824a
# ╠═97d9b1b4-453b-4e2b-a7dd-4dc941900c00
# ╠═c13c7337-db5d-4177-8a6c-79ff3c1173f2
# ╠═db80c34b-7970-4d72-89b4-c90091c1f553
# ╠═03bae442-29ee-4272-b630-d051956a60ab
# ╟─869ec51c-e6f5-4733-84f8-6b006dc0f294
# ╠═8d99a321-f9c6-4da9-921d-55bf50034dc5
# ╠═01aae5a0-f208-477d-974b-3945f0d56fd0
# ╠═ab81fd7e-f6fe-4c25-b5b1-9f19f4d9cddb
# ╠═245afe84-62df-4496-b8ff-c4d586fb0f47
# ╠═db717a7f-bc68-47fb-b894-ff97bd1bd1ec
# ╟─92248a65-e3fe-4071-b410-ea853b4c4c11
# ╠═bc322aa5-c57f-487c-9e5a-c907e3582ef7
# ╠═1621bce6-e8bb-4d70-8784-21aa80474af3
# ╠═1dc3eec5-fe44-474a-a8df-007102a522d3
# ╠═dadb9d11-80ba-4e85-b49d-e7671672dc6e
# ╠═983d0f16-1484-44c1-818d-5827ab582e5f
# ╠═5abe2ecd-51bf-4854-b287-1ce85ab30347
# ╠═192c4600-a9c1-4dbb-a74c-6352bb125059
# ╠═81be70da-010d-45b6-9e8e-d5500dcbb867
# ╠═c80a9688-3ceb-45e2-8c18-93d0599f03fc
# ╠═bbdecd68-45cc-4894-8d19-cd13a51bdedc
# ╠═bd33b46b-bfdc-449e-a614-cc04688f5b89
# ╠═58e67962-6cec-4d01-b4f8-ce003c0e4967
# ╠═9cdcd137-a2fb-4f82-98a3-e613e1a80f54
# ╠═c96c018b-d4a4-4ebf-925c-bbf9bcb4fc8f
# ╠═f6b3bb75-3565-41f0-b592-81b2b36eb651
# ╟─35f82e76-05c0-4aa4-823e-16187d90ab25
# ╟─9de81723-c470-4e53-8fb8-ea2c1a889c1f
# ╠═2b1689f4-89a7-42a3-916d-5731dc91b5b3
# ╠═d5f7057c-da96-405e-9ad8-a8032c7c47a1
# ╠═926897bd-edca-4961-975b-7117b7d000e9
# ╠═07569676-36be-4b54-8105-c28e97bb872a
# ╠═122767a3-3348-4b4c-990f-315dfe7c30ca
# ╠═156e04fd-cf99-4e11-97a9-eb334a85a2cb
# ╠═636135ce-bfd0-4786-ad99-25af1832fe2f
# ╠═ea7dd844-8d79-4d47-98ce-c91fd3fa04b6
# ╟─55317d1e-88de-4b61-8b1a-616463ab7299
# ╟─a4a7ecb8-293b-47a4-a45a-b0fbc7f9fc7e
# ╠═0bf2dc1c-084f-4714-8bb6-11e631767199
# ╠═213e99b1-13ba-48cb-9838-f80077d85992
# ╠═d2773d45-b995-48b4-bcf8-754341dc3dcd
# ╟─ad624619-cc2e-42c5-ae62-817771f600ee
# ╠═951cde36-9000-4d8f-9cb5-81401be51671
# ╟─94daad38-215e-4f54-8165-af5f8a373804
# ╟─aa3c015d-371a-444d-9b9c-98d2aef53d1e
# ╠═990a8734-a7a9-4d11-b10f-447912e47e60
# ╠═fe573923-fd00-4c76-ac58-eedac559b6a6
# ╠═38038448-d16d-4e0c-8ac1-233aca2635e3
# ╠═c0a4b035-baf1-452f-a608-95c4eeaa6549
# ╠═d7666993-dc60-447c-ab90-696b2c9862c1
# ╠═649d0694-704f-4de5-a882-338a276968c8
# ╠═9fabd4e0-b733-48e2-81bf-b0c6b5a4dae1
# ╠═1945ba42-eb1d-469f-bf53-c5b8e7573b54
# ╟─4a4d0612-779f-4715-beb6-913f8023cccd
# ╠═5f46bdeb-31af-417e-9e80-291898252734
# ╠═7c2560f8-c8cc-4658-ab17-77d8322b2dfd
# ╟─a2b2a498-524c-47f2-a15f-bdc0768f820a
# ╠═95eaff97-3e79-4b3b-8265-0bd77bd31183
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
