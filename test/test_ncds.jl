function test_ncds(file_name)
    inst= Instance(file_name,0)
    solgroup=OT_group(inst)
    solgrouprelax=OT_group(inst,0.2, 0.199)
    soljoint=OT_joint(inst, 0.0, 0.0, 0.2)
    soljointrelax=OT_joint(inst, 0.199, 0.1, 0.2)
    effjoint = round.(Int,inst.nB*soljoint.jointYZB+inst.nA*soljoint.jointYZA);
    effjointrelax = round.(Int,inst.nB*soljointrelax.jointYZB+inst.nA*soljointrelax.jointYZA);
    effgroup = round.(Int,inst.nB*solgroup.jointYZB+inst.nA*solgroup.jointYZA);
    effgrouprelax = round.(Int,inst.nB*solgrouprelax.jointYZB+inst.nA*solgrouprelax.jointYZA);
    tab = [sum((inst.Yobserv.== y) .& (inst.Zobserv.==z)) for y in inst.Y, z in inst.Z];
    println("Error with group = ", 1/2.0*sum(abs.(tab.-effgroup))/sum(tab));
    println("Error with relaxed group = ",1/2.0*sum(abs.(tab.-effgrouprelax))/sum(tab));
    println("Error with  joint = ",1/2.0*sum(abs.(tab.-effjoint))/sum(tab));
    println("Error with relaxed joint = ",1/2.0*sum(abs.(tab.-effjointrelax))/sum(tab));
    return tab,effgroup,effgrouprelax,effjoint,effjointrelax
end
