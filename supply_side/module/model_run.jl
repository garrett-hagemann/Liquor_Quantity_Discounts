# addprocs(68) # for Stampede 2 on 1 node
addprocs(2) # for home computer testing

@everywhere srand(69510606) #seeding random number gen

@everywhere include("LiquorQD.jl")
@everywhere include("DataSetup.jl")
@everywhere using DataSetup, LiquorQD

markets_array = data_setup()

@everywhere function mkt_est(m::Market)
    out = Dict{Tuple{Market,Liquor},Array{Float64,2}}()
    println("Working with ", m.year,"-",m.month,".")
    for j in m.products
        println("working with product ", j.id)
        if !isnull(j.ps)
            println("Matched price schedule.")
          tmp_ps = get(j.ps) # becase the ps field is nullable, need to use get
          print("Generating perturbed price schedules. ")
          dev_ps = dev_gen(tmp_ps,0.05)
          println("Done.")
          print("Pre-calculating retail prices. ")
          pre_calc = Dict{Int64,Float64}[]
          for s in dev_ps
            tmp_dict = Dict(i => p_star(s.rhos[i],j,nested_coefs,obs_inc_dist,m) for i in 1:(s.N-1))
            push!(pre_calc,tmp_dict)
          end
          println("Done.")
          print("Pre-calculating shares at retail prices. ")
          s_pre_calc = Dict{Int64,Float64}[]
          for d in pre_calc
              tmp_s_dict = Dict(key=>share(p,j,nested_coefs,obs_inc_dist,m) for (key,p) in d)
              push!(s_pre_calc,tmp_s_dict)
          end
          println("Done.")
          min_rho = minimum(tmp_ps.rhos)
          sol,xtrace,ftrace = optimize_moment(tmp_ps,dev_ps,j,nested_coefs,obs_inc_dist,m,25,pre_calc,s_pre_calc,x0=[min_rho/2.0,1.0])
          println(sol)
          trace = [vcat(xtrace'...) ftrace]
          out[(m,j)] = trace
        else
            println("No matching price schedule data.")
        end
    end
    return out
end

mkts_for_est = markets_array
res = pmap(mkt_est,mkts_for_est)

out_dict = merge(res...)

for (key,v) in out_dict
    m = key[1]
    j = key[2]
    out_str = "traces/trace_$(m.year)_$(m.month)_$(j.id).csv"
    writedlm(out_str,v)
end



#=
m = markets_array[1]
j = m.products[1]
tmp_ps = get(j.ps)
println("Generating deviation price schedules.")
dev_ps = dev_gen(tmp_ps,0.2)
println("Pre-calculating retail prices.")
pre_calc = Dict{Int64,Float64}[]
for s in dev_ps
  tmp_dict = Dict(i => p_star(s.rhos[i],j,coef_array,inc_weights,m) for i in 1:(s.N-1))
  push!(pre_calc,tmp_dict)
end

println("b, Q")
for b = 0.5:.5:20
  tmp_params = WholesaleParams(8.0,1.0,b)
  res = moment_obj_func(tmp_ps,dev_ps,tmp_params,j,coef_array,inc_weights,m,pre_calc)
  println(b,",",res,"\n")
end
=#
