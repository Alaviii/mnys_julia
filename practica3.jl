using DelimitedFiles

function main()
    sistemagauss::Array{Float64, 2} = readdlm("matrizin.dat")
    sistemadiv::Array{Float64, 2} = readdlm("matrizin.dat")
    sistemajordan::Array{Float64, 2} = readdlm("matrizin.dat")
    matrizcuadrada::Array{Float64, 2} = readdlm("matrizcuadrada.dat")
    sistema::Array{Float64,2} = readdlm("matrizin.dat")
    solvible::Array{Float64,2} = readdlm("sistemasolvible.dat")


    writedlm("matrizout.dat",sistemagauss)

    display(gauss!(sistemagauss))

    display(divisionunica!(sistemadiv))

    display(gaussjordan!(sistemajordan))

    (l,u) = lu(matrizcuadrada)
    display(l)
    display(u)
    display(l*u)

    display(richardson(sistema,iteraciones = 150))
    display(jacobi(sistema,iteraciones = 150))
    display(gaussseidel(sistema,iteraciones = 150))

    display(richardson(solvible,iteraciones = 1500))
    display(jacobi(solvible,iteraciones = 1500))
    display(gaussseidel(solvible, errorlimite = 0.01))

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

function lu(a::Array{Float64,2})
    if size(a,1) != size(a,2) return "No es una matriz cuadrada", "No es una matriz cuadrada"

    else
        n = size(a,1)
        l = zeros(Float64,n,n)
        u = zeros(Float64,n,n)
        for i in range(1,n)
            l[i,i] = 1
        end
        
        for p in range(1,n)
            for j in range(p,n)
                suma = 0

                for k in range(1,p-1)
                    suma += l[p,k]*u[k,j]
                end

                u[p,j] = a[p,j] - suma
            end

            for i in range(p+1,n)
                suma = 0

                for k in range(1,p-1)
                    suma += l[i,k]*u[k,p]
                end

                l[i,p] = (a[i,p]-suma)/u[p,p]
            end
        end

        return l,u
    end
end

function richardson(a::Array{Float64,2}; seed::Array{Float64,1}=zeros(Float64,size(a,1)), iteraciones = 100000, errorlimite = 0)
    if size(a,1)+1 != size(a,2) return "No es una matriz ampliada"

    else
        n = size(a,1)
        indice_b = size(a,2)

        xvieja = seed
        xnueva = zeros(Float64,size(a,1))

        k = 0

        error = errorlimite+1

        while (k<iteraciones && error>errorlimite)
            for i in range(1,n)
                suma = 0
                for j in range(1,n)
                    suma += a[i,j]*xvieja[j]
                end
                xnueva[i] = xvieja[i] - suma + a[i,indice_b]
            end

            k += 1

            error = 0

            for i in range(1,n)
                if abs(xnueva[i]-xvieja[i])>error
                    error = abs(xnueva[i] - xvieja[i])
                end
            end

            xvieja = copy(xnueva)
        end
        
        return xnueva
    end
end

function jacobi(a::Array{Float64,2}; seed::Array{Float64,1}=zeros(Float64,size(a,1)), iteraciones = 100000, errorlimite = 0)
    if size(a,1)+1 != size(a,2) return "No es una matriz ampliada"

    else
        n = size(a,1)
        indice_b = size(a,2)

        xvieja = seed
        xnueva = zeros(Float64,size(a,1))

        k = 0

        error = errorlimite+1

        while (k<iteraciones && error>errorlimite)
            for i in range(1,n)
                suma = 0
                for j in range(1,n)
                    if i!=j
                        suma += a[i,j]*xvieja[j]
                    end
                end

                xnueva[i] = (a[i,indice_b]-suma)/a[i,i]
            end

            k += 1

            error = 0

            for i in range(1,n)
                if abs(xnueva[i]-xvieja[i])>error
                    error = abs(xnueva[i] - xvieja[i])
                end
            end

            xvieja = copy(xnueva)
        end
        
        return xnueva
    end
end

function gaussseidel(a::Array{Float64,2}; seed::Array{Float64,1}=zeros(Float64,size(a,1)), iteraciones = 100000, errorlimite = 0)
    if size(a,1)+1 != size(a,2) return "No es una matriz ampliada"

    else
        n = size(a,1)
        indice_b = size(a,2)

        xvieja = seed
        xnueva = zeros(Float64,size(a,1))

        k = 0

        error = errorlimite+1

        while (k<iteraciones && error>errorlimite)
            for i in range(1,n)
                sumanueva = 0
                for j in range(1,i-1)
                    sumanueva += a[i,j]*xnueva[j]
                end
                sumavieja = 0
                for j in range(i+1,n)
                    if j != i
                        sumavieja += a[i,j]xvieja[j]
                    end
                end

                xnueva[i] = (a[i,indice_b]-sumanueva-sumavieja)/a[i,i]
            end
            
            k += 1

            error = 0

            for i in range(1,n)
                if abs(xnueva[i]-xvieja[i])>error
                    error = abs(xnueva[i] - xvieja[i])
                end
            end

            xvieja = copy(xnueva)
        end

        return xnueva
    end
end


main()