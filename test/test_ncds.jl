@testset "NCDS" begin 

    using OTRecod

    inst          = Instance(file_name, 0)
    solgroup      = ot_group(inst)
    solgrouprelax = ot_group(inst,0.2, 0.199)
    soljoint      = ot_joint(inst, 0.0, 0.0, 0.2)
    soljointrelax = ot_joint(inst, 0.199, 0.1, 0.2)
    effjoint      = round.(inst.nB * soljoint.jointYZB      + inst.nA * soljoint.jointYZA)
    effjointrelax = round.(inst.nB * soljointrelax.jointYZB + inst.nA * soljointrelax.jointYZA)
    effgroup      = round.(inst.nB * solgroup.jointYZB      + inst.nA * solgroup.jointYZA)
    effgrouprelax = round.(inst.nB * solgrouprelax.jointYZB + inst.nA * solgrouprelax.jointYZA)

    tab = [sum((inst.Yobserv.== y) .& (inst.Zobserv.==z)) for y in inst.Y, z in inst.Z]

    println(" Error with group         = ", 1/2.0 * sum(abs.(tab.-effgroup))/sum(tab))
    println(" Error with relaxed group = ", 1/2.0 * sum(abs.(tab.-effgrouprelax))/sum(tab))
    println(" Error with  joint        = ", 1/2.0 * sum(abs.(tab.-effjoint))/sum(tab))
    println(" Error with relaxed joint = ", 1/2.0 * sum(abs.(tab.-effjointrelax))/sum(tab))

end
