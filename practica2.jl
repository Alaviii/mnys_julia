using DelimitedFiles

cotes::Array{Int64,2} = readdlm("cotes.dat")

function main()
    function f(x)
        return (sqrt(8*0.0077))/sqrt(2.3^2.5-x^2.5)
    end

    function fbeta(x, beta)
        return sqrt(8*0.0077)/sqrt(2.3^beta-x^beta)
    end

    limsup = 2.3-0.001

    display(newtoncotes(f,0,limsup,6))
    display(trapeciocompuesto(f,0,limsup,20))
    display(simpsoncompuesto(f,0,limsup,20))
    display(romberg(f,0,limsup,6))

    display(segundaderivada(f,1.77,0.001))

    beta = 2
    n = (3-2)/100
    suma = 0
    i = 0
    while beta<=3
        suma += trapeciocompuestobeta(fbeta,beta,0,limsup,20)
        beta+=n
        i +=1
    end

    display(suma)
    display(i)

end

function newtoncotes(func::Function, lim_inf, lim_sup, n_newton_cotes::Int64)
    integral::Float64 = 0

    h::Float64 = (lim_sup-lim_inf)/n_newton_cotes

    for i in UnitRange(1,n_newton_cotes+1)
        integral += (cotes[n_newton_cotes,i]/cotes[n_newton_cotes,10])*func(lim_inf+(i-1)*h)
    end
    
    return (lim_sup-lim_inf)*integral
end

function trapeciocompuesto(func::Function, lim_inf, lim_sup, n_intervalos::Int64)
    h::Float64 = (lim_sup-lim_inf)/n_intervalos

    integral::Float64 = (h/2)*(func(lim_inf)+func(lim_sup))

    suma::Float64 = 0

    for i in UnitRange(1,n_intervalos-1)
        suma += func(lim_inf + i*h)
    end
    
    return integral + h*suma
end

function trapeciocompuestobeta(func::Function, beta, lim_inf, lim_sup, n_intervalos::Int64)
    h::Float64 = (lim_sup-lim_inf)/n_intervalos

    integral::Float64 = (h/2)*(func(lim_inf,beta)+func(lim_sup,beta))

    suma::Float64 = 0

    for i in UnitRange(1,n_intervalos-1)
        suma += func(lim_inf + i*h,beta)
    end
    
    return integral + h*suma
end

function simpsoncompuesto(func::Function, lim_inf, lim_sup, n_intervalos::Int64)
    if n_intervalos%2 != 0
        return "error, n ha de ser par"
    end
    
    h::Float64 = (lim_sup-lim_inf)/n_intervalos

    suma_pares::Float64 = 0
    suma_impares::Float64 = 0

    i::Int64 = 1

    while i<n_intervalos
        suma_impares += func(lim_inf + i*h)
        i+= 2
    end

    i = 2

    while i<n_intervalos
        suma_pares += func(lim_inf + i*h)
        i+= 2
    end
    
    return (h/3)*(func(lim_inf)+func(lim_sup)+4*suma_impares+2*suma_pares)
end

function romberg(func::Function, a, b, p::Int64)
    if p == 0
        return ((b-a)/2)*(func(a)+func(b))
    else
        suma = 0
        h2p = (b-a)/(2^p)

        for i in range(1,2^(p-1))
            suma += func(a+(2*i-1)*h2p)
        end

        return 0.5*romberg(func,a,b,p-1) + h2p*suma
    end
end

function derivadadefinicion(func::Function, x, h)
    return (func(x+h)-func(x))/h
end

function derivadataylor(func::Function, x, h)
    return (func(x+h)-func(x-h))/(2*h)
end

function derivadarichardson(f::Function, x, h)
    return (4.0/3.0)*((f(x+h/2)-f(x-h/2))/h) - (1.0/3.0)*((f(x+h)-f(x-h))/(2*h))
end

function segundaderivada(f::Function, x, h)
    return (f(x+h)-2*f(x)+f(x-h))/h^2
end

main()
