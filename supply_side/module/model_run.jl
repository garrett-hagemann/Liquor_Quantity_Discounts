# addprocs(68) # for Stampede 2 on 1 node
#addprocs(2) # for home computer testing

data_file_arg = ARGS[1] # needs a command line argument passed. Is the data file

@everywhere srand(69510606) #seeding random number gen

@everywhere include("LiquorQD.jl")
@everywhere include("DataSetup.jl")
@everywhere using DataSetup, LiquorQD

markets_array = data_setup(data_file_arg)

@everywhere function mkt_est(tup::Tuple{Market,Liquor})
    m = tup[1]
    j = tup[2]
    out = Dict{Tuple{Market,Liquor},Array{Float64,2}}() # out should be a 1 element Dict. Want Dict because we merge Dicts later
    println("working with product $(j.id) in $(m.year)-$(m.month)")
    if !isnull(j.ps)
        println("Matched price schedule.")
      tmp_ps = get(j.ps) # becase the ps field is nullable, need to use get
      # hack-y code
      println(tmp_ps.t_cuts)
      scale = (1-1/tmp_ps.N)/tmp_ps.t_cuts[end]
      tmp_ps.t_cuts = tmp_ps.t_cuts.*scale
      println(tmp_ps.t_cuts)
      #end hack
      print("Generating perturbed price schedules. ")
      dev_ps = dev_gen(tmp_ps,0.025)
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
      sol,xtrace,ftrace = optimize_moment(tmp_ps,dev_ps,j,nested_coefs,obs_inc_dist,m,500,pre_calc,s_pre_calc,x0=[min_rho/2.0,1.0])
      #println(sol)
      trace = [vcat(xtrace'...) ftrace]
      #trace = [sol.c sol.b ftrace]
      out[(m,j)] = trace
    else
        println("No matching price schedule data.")
    end
    return out
end

mkts_for_est = [(m,j) for m in markets_array for j in m.products] # array of tuples
res = pmap(mkt_est,mkts_for_est)

for d in res
    if !isempty(d)
        for (key,v) in d
            m = key[1]
            j = key[2]
            out_str = "traces/trace_$(m.year)_$(m.month)_$(j.id).csv"
            writecsv(out_str,v)
        end
    end
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
