
""" Tests on the code """

include("main.jl")

function test(instance_path::String)

    println("-------------------------------------------")
    println("INSTANCE : ",instance_path)
    println(CONFIG)

    dichos = [JUMP, COMBO]
    methods = [EXACT_JUMP, EXACT_COMBO, RELAX_LIN_CLASSIC, RELAX_LIN_SPEED_UP]
    tab = [false, true]

    last_cpt, last_lb = 0,0

    for dicho in dichos
        COMPONENTS.firstDicho = dicho
        for method in methods
            COMPONENTS.methodUB = method
            for vN in tab
                COMPONENTS.nadirsShift = vN
                for vH in tab
                    COMPONENTS.primal_heuristic = vH
                    cpt, lb = main(instance_path)
                    #cpt, lb = main()
                    if last_cpt != 0
                        @assert lb.sols == last_lb.sols "Different results"
                        println("   SUCCESS - ",COMPONENTS)
                    else
                        last_cpt = cpt
                        last_lb = lb
                        println("   FIRST   - ",COMPONENTS)
                    end
                end
            end
        end
    end
end

function test()
    for size in [10, 20, 30, 40, 50]
        for id in 1:5
            path = string("instances/Inst_",size,"_",id,".dat")
            test(path)
        end
    end
end