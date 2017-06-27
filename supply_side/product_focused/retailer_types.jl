#= This file recovers the types for each (p,q) available in the dataset.
The approach is as follows:
1) Figure out the cutoff types associated with the cutoff quantities in the price
schedule
2) given the cutoff types, calculate the prices the cutoff types would charge
3) Retailers who charge prices between those will be in that segment of the price schedule
4) Given the price segment the retailer faces, can back out the retailer type.

Cavets:
1) using the demand estimates from the berry logit. So no random coefs. This will be more complicated
in the presence of random coefficients
2) This uses the market size as defined by the number of trips in the nielsen data. This can also
be adjusted, but the demand estimation would need to be adjusted as well
3) This relies on the assumption that the optimal price is monotone in the retailer type across price segments.
It is clearly monotone within a price segment, but across segments might not be. However, note that this would
imply that the optimal price function is discontinuous in type. That might be able to be checked.=#

# Packages to include
using DataArrays, DataFrames, Roots, ReverseDiffSource

df = readtable("../../demand_estimation/berry_logit/berry_logit.csv")
char_list = [:price, :d_gin, :d_vod, :d_rum, :d_sch, :d_brb, :d_whs, :holiday]
coefs = readtable("../../demand_estimation/berry_logit/blogit_coefs.csv")

#= each row in df represents one bundle of p,q. Each row also contains all
the needed information to determine the cutoffs. So, we'll be looping over the relevant
price schedule segments and cutoffs and doing it for each row. =#

#= each column in coefs is the appropriate coefficient =#
coef_vec = convert(Array,coefs) # row vector
alpha = coef_vec[1,1]
#= We're going to just define the share function, then use the automatic differentiation
package to get the required derivatives to get the correct equations for the type and optimal
price. This saves us redefining things repeatedly and doing derivatives by hand =#

# try for prod = 650, mkt = 178. Ultimately want to loop over both
# getting product characteristics

markets = convert(Vector, levels(df[:,:mkt])) # needs to be a vector to play nice with rdiff

csvfile = open("retailer_types.csv", "w")
write(csvfile, "product,mkt,lambda\n")

for market in markets
	products = levels(df[df[:mkt] .== market,:product])
	
	for product in products
		println("Working with Market: $market, Product: $product")

		# defining selection boolians so that we can select appropriate products/markets
		prod_bool = ((df[:product] .== product) & (df[:mkt] .== market))
		other_bool = ((df[:product] .!= product) & (df[:mkt] .== market))

		prod_chars = convert(Array,df[prod_bool,char_list])
		other_chars = convert(Array,df[other_bool,char_list])

		# flag that determines if product has matched price data
		matched_upc = (df[prod_bool, :_merge_purchases][1] == 3)

		#Defining some constants for use in the share function
		global xb_prod = (prod_chars[:,2:end]*coef_vec[:,2:end]')[1] + df[prod_bool,:xi][1] # should be scalar
		if size(other_chars)[1] == 0
			global xb_others = 0
		else
			global xb_others = other_chars*coef_vec' + convert(Vector,df[other_bool,:xi]) # should be J-1 x 1
		end

		#Grabbing market size. Needs to be the same as that used to estimate shares
		M = df[prod_bool,:M][1]
		prod_price = df[prod_bool, :price][1]


		#Defining share function which is ONLY a function of price
		function share(p)
			#= Function calculates shares of a given  product in a given market 
			based on price. Characteristics of the product and all other products
			are assumed fixed. 

			p = price =#
			
			num = exp(alpha*p + xb_prod)
			denom = 1 + num + sum(exp(xb_others))
			
			share = num/denom	

			return share
		end

		#= Derivative of the share fuction. Returns a tuple that is (f(x), f'(x)) =#
		d_share = rdiff(share, (10,), order = 1) # the 10 is arbitrary. Just needed for type of function argument

		#= For products that have matched price schedule data we want to figure out the retail price 
		that will result in the cutoff quantity being sold for each part of the price schedule. That is,
		for a price schedule with cutoffs at 5 and 10 cases, we want to find a p such that s(p)*M = 5*12 (with
		12 bottles per case) and p such that s(p)*M = 10*12. This can be done easily by finding the root of 
		the function g(p) = s(p) - cutoff. 

		With this information, we can determine what part of the price schedule the retailer
		is facing (because their price must be between the cutoff prices for the next tier. 
		This allows us to figure out the retailer type, the ultimate goal =#


		if matched_upc == true
			println("Product has matched price data. Evaluating cutoff prices.")
			
			if df[prod_bool,:discount_size_type][1] == "BOTTLE"
				bot_per_case = 1
			else
				bot_per_case = df[prod_bool,:bottles_per_case][1]
			end
			
			cutoff_q_list = [:disc_q0, :disc_q1, :disc_q2, :disc_q3, :disc_q4, :disc_q5, :disc_q6,
				:disc_q7, :disc_q8, :disc_q9, :disc_q10]
			
			cutoff_p_list = []
			
			rho_list = [:actual_p0, :actual_p1, :actual_p2, :actual_p3, :actual_p4, :actual_p5, :actual_p6,
				:actual_p7, :actual_p8, :actual_p9, :actual_p10]

			for var in cutoff_q_list
				cutoff_q = df[prod_bool, var][1]
				if isna(cutoff_q) # indicates that product has no discount info for this cut
					break
				else
					g(p) = d_share(p)[1]*M - cutoff_q*bot_per_case #d_share[1] is the value of the function. NOT SURE why this works here and not regular share function?
					cutoff_p = fzero(g,0,1e9) #uses derivative free method
					push!(cutoff_p_list,cutoff_p)
				end
			end
			bool_p = (prod_price .> cutoff_p_list) # boolean list
			seg_ind = find(bool_p) # finds the first True element. Returns 0 if list is all false
			if length(seg_ind) == 0
				rho = df[prod_bool,rho_list[length(cutoff_p_list)]][1]
			elseif seg_ind[1] == 1
				rho = df[prod_bool,rho_list[1]][1]
			else
				rho = df[prod_bool,rho_list[seg_ind[1]-1]][1]
			end

			#= With the part of the price schedule this retailer is pricing from,
			we can now determine their type =#
			
			lambda = prod_price - rho + share(prod_price)/d_share(prod_price)[2] # 2 index is the derivate of share function
			out_tuple = (product, market, lambda)
			write(csvfile, join(out_tuple,","),"\n")
		else
			println("Product has no matching price data.")
		end
	end
end
