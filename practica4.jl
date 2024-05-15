using DelimitedFiles

function main()
    function exmasx(x)
        return ℯ^x + x
    end
    function polinomiomalvado(x)
        return x^3-5x^2+7x-3
    end

    display(biseccion(exmasx, -1, 0, 10^-3))

    display(secante(exmasx,0,1,10^-3))
    
end

function biseccion(f::Function,liminf,limsup,errorlimite)
    if f(liminf)*f(limsup) > 0
        return "No se cumple el teorema de bolzano, otro intervalo"
    else
        a = liminf
        b = limsup

        if f(a) == 0 
            return a
        elseif f(b) == 0 
            return b
        end

        iteration = 0
        c = 0 #no se por qué pero peta si no
        
        error = f((a+b)/2)

        while abs(error)>abs(errorlimite)
            c = (a+b)/2

            if abs(f(a)*f(c)) < 0.00000000000001
                return c
            elseif f(a)*f(c) < 0
                b = c
            else 
                a = c
            end

            iteration += 1

            error = f(c)
        end

        return c
    end
end

function falsaposicion(f::Function,liminf,limsup,errorlimite)
    if f(liminf)*f(limsup) > 0
        return "No se cumple el teorema de bolzano, otro intervalo"
    else
        a = liminf
        b = limsup

        if f(a) == 0 
            return a
        elseif f(b) == 0 
            return b
        end

        iteration = 0
        c = 0
    end
end

function newtonraphson(f::Function,fprima::Function,semilla,errorlimite)
    xnow = semilla

    iteration = 0

    error = errorlimite+1

    xpost = 0 #npi

    while abs(error)>abs(errorlimite)
        xpost = xnow - f(xnow)/fprima(xnow)

        error = f(xpost)

        xnow = xpost

        iteration += 1
    end

    return xpost
end

function secante(f::Function,semilla0,semilla1,errorlimite)
    xpre = semilla0
    xnow = semilla1

    iteration = 0

    error = errorlimite+1

    xpost = 0 #npi

    while abs(error)>abs(errorlimite)
        xpost = xnow - f(xnow)*((xnow-xpre)/(f(xnow)-f(xpre)))

        error = f(xpost)

        xpre = xnow
        xnow = xpost

        iteration += 1
    end

    return xpost
end



main()