using DelimitedFiles

function minimoscuadrados(tabla::Array{Float64,2})
    N::Int64 = size(tabla,1)-1

    sumax::Float64 = sumay::Float64 = sumaxy::Float64 = sumaxx::Float64 = 0

    for i in UnitRange(1,N+1)
        sumax += tabla[i,1]
        sumay += tabla[i,2]
        sumaxy += tabla[i,1]*tabla[i,2]
        sumaxx += tabla[i,1]*tabla[i,1]
    end

    a1::Float64 = (((N+1)*sumaxy-sumax*sumay)/((N+1)*sumaxx-sumax*sumax))

    a0::Float64 = (sumaxx*sumay-sumaxy*sumax)/((N+1)*sumaxx-sumax*sumax)

    return (a0, a1)
end

function lagrange(tabla::Array{Float64,2}, valor::Float64)
    n::Int64 = (size(tabla,1)-1)
    lagrange_result::Float64 = 0

    for i in UnitRange(1,n+1)
        li = 1
        for j in UnitRange(1,n+1)
            if(j!=i)
                top = valor - tabla[j,1]
                bot = tabla[i,1] - tabla[j,1]
                li *= top/bot
            end
        end
        lagrange_result += tabla[i,2]*li
    end

    return lagrange_result
end



function main()
    tabla::Array{Float64,2} = readdlm("tabla.dat")

    fichout::IOStream = open("test.dat", "w")

    lagrange_input = 0.4
    while lagrange_input <=0.8
        println(fichout, lagrange(tabla,lagrange_input))

        lagrange_input +=0.01
    end
end

main()