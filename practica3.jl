using DelimitedFiles

function main()
    sistemagauss::Array{Float64, 2} = readdlm("matrizin.dat")
    sistemadiv::Array{Float64, 2} = readdlm("matrizin.dat")
    sistemajordan::Array{Float64, 2} = readdlm("matrizin.dat")
    matriz::Array{Float64, 2} = readdlm("matrizcuadrada.dat")


    writedlm("matrizout.dat",sistemagauss)

    display(gauss!(sistemagauss))
    display(divisionunica!(sistemadiv))
    display(gaussjordan!(sistemajordan))
    display(lu(matriz))
end

function gauss!(sistema::Array{Float64,2})
    if size(sistema,1)+1 != size(sistema,2) return "No es una matriz ampliada"

    else
        n = size(sistema,1)
        for iteration in range(1,n-1)
            pivot = sistema[iteration,:]
            for selectedrow in range(iteration+1,n)
                sistema[selectedrow,:] -= pivot*(sistema[selectedrow,iteration]/pivot[iteration])
            end
        end

        coeficientes = Float64[]
        resize!(coeficientes,n)

        for i in range(start=n,stop=1,step=-1)
            suma = 0
            for j in range(i+1,n)
                suma += sistema[i,j]*coeficientes[j]
            end
            coeficientes[i] = ((sistema[i,n+1] - suma))/ sistema[i,i] 

        end

        return coeficientes
    end
end

function divisionunica!(sistema::Array{Float64,2})
    if size(sistema,1)+1 != size(sistema,2) return "No es una matriz ampliada"

    else
        n = size(sistema,1)
        for iteration in range(1,n)
            sistema[iteration,:] /= sistema[iteration,iteration]
            pivot = sistema[iteration,:]
            for selectedrow in range(iteration+1,n)
                sistema[selectedrow,:] -= pivot*sistema[selectedrow,iteration]
            end
        end

        coeficientes = Float64[]
        resize!(coeficientes,n)

        for i in range(start=n,stop=1,step=-1)
            suma = 0
            for j in range(i+1,n)
                suma += sistema[i,j]*coeficientes[j]
            end
            coeficientes[i] = ((sistema[i,n+1] - suma))

        end

        return coeficientes
    end
end

function gaussjordan!(sistema::Array{Float64,2})
    if size(sistema,1)+1 != size(sistema,2) return "No es una matriz ampliada"

    else
        n = size(sistema,1)
        for iteration in range(1,n)
            sistema[iteration,:] /= sistema[iteration,iteration]
            pivot = sistema[iteration,:]
            for selectedrow in range(1,n)
                if selectedrow != iteration
                    sistema[selectedrow,:] -= pivot*sistema[selectedrow,iteration]
                end
            end
        end

        coeficientes = Float64[]
        resize!(coeficientes,n)

        for i in range(start=n,stop=1,step=-1)
            coeficientes[i] = sistema[i,n+1]
        end

        return coeficientes
    end
end

function lu(matriz::Array{Float64,2})
    if size(matriz,1) != size(matriz,2) return "No es una matriz ampliada"

    else
        n = size(matriz,1)

        l = zeros(Float64,n,n)
        u = zeros(Float64,n,n)

        for i in range(1,n)
            l[i,i] = 1
        end

        return l,u
    end
end

main()