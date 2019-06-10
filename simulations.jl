dirname("C:/Users/Kiki/Google Drive/3_Recherche/R1_Recherche/5_INSA/DataMerging_4/Programmes/mainDir/K-0.5/")
dirname("C:/Users/vagares/Google Drive/3_Recherche/R1_Recherche/5_INSA/DataMerging_4/Programmes/mainDir/K-0.5/")
using JuMP
using Gurobi
using Ipopt
using Base

M=1000
simul=[]
results= spzeros(M,nf) 
for i in 1:M :
  data_file="C:/Users/Kiki/Google Drive/3_Recherche/R1_Recherche/5_INSA/DataMerging_4/Programmes/mainDir/K-0.5/tab$(i).txt"
  data_file="C:/Users/vagares/Google Drive/3_Recherche/R1_Recherche/5_INSA/DataMerging_4/Programmes/mainDir/K-0.5/tab$(i).txt"

  inst = Instance(data_file)
  
  # solve models
  # individual transports
  sol_indiv_base = solve_indiv_base(inst, false, true)
  sol_indiv_relax = solve_indiv_relax(inst)
  sol_indiv_reg = solve_indiv_reg(inst, norm2, 2, 50.0)
  
  # group transport
  sol_group_base = solve_group_base(inst,optimal)
  sol_group_relax = solve_group_relax(inst, 0.5,optimal)
  
  estim_indiv_base = compute_pred_error(inst,sol_indiv_base.YAtrans,sol_group_base.YBtrans)
  estim_indiv_relax = compute_pred_error(inst,sol_indiv_relax.YAtrans,sol_group_base.YBtrans)
  estim_indiv_reg = compute_pred_error(inst,sol_indiv_reg.YAtrans,sol_group_base.YBtrans)
  estim_group_base = compute_pred_error(inst,sol_group_base.YAtrans,sol_group_base.YBtrans)
  estim_group_relax = compute_pred_error(inst,sol_group_relax.YAtrans,sol_group_base.YBtrans)
  
  simul= [simul; estim_indiv_base estim_indiv_relax  estim_indiv_reg estim_group_base estim_group_relax]
end

