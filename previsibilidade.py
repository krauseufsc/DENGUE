from math import *

def soma_lista(lista1, lista2): #soma os valores paralelos da lista e retorna essa soma em uma lista de mesmo tamanho
    saida = []
    for i in range(len(lista1)):
        soma = lista1[i] + lista2[i]
        saida.append(soma)
    return saida

def multiplica_lista(const, lista): #multiplica uma constante a todos os valores de uma lista
    saida = []
    for i in lista:
        aux = i*const
        saida.append(aux)
    return saida

def nega_zero(lista): #transforma todos os valores negativos de uma lista em zero
    nova_lista = []
    for i in lista:
        if i < 0:
            i = 0
        nova_lista.append(i)
    return nova_lista

def field(inp,t,n):
    s = inp[0]
    i1 = inp[1]
    i2 = inp[2]
    r1 = inp[3]
    r2 = inp[4]
    s1 = inp[5]
    s2 = inp[6]
    i12 = inp[7]
    i21 = inp[8]
    r = inp[9]
    sv = inp[10]
    v1 = inp[11]
    v2 = inp[12]

    m = sv+v1+v2
    phi=0.8
    mu=1/65
    alfa=2
    gama=52
    nu=36.5
    teta=2*nu
    omega=2*3.14*6
    xi=nu*(1+0.4 * cos(omega*t))
    beta=2*gama

    ds = -(beta/m)*s*(v1+v2)+mu*(n-s)
    di1 = (beta/m)*s*v1-(gama+mu)*i1
    di2 = (beta/m)*s*v2-(gama+mu)*i2
    dr1 = gama*i1-(alfa+mu)*r1
    dr2 = gama*i2-(alfa+mu)*r2
    ds1 = -(beta/m)*s1*v2+alfa*r1-mu*s1
    ds2 = -(beta/m)*s2*v1+alfa*r2-mu*s2
    di12 = (beta/m)*s1*v2-(gama+mu)*i12
    di21 = (beta/m)*s2*v1-(gama+mu)*i21
    dr = gama*(i12+i21)-mu*r
    dsv = -(teta/n)*sv*(i1+i2+phi*(i12+i21))+xi*m-nu*sv
    dv1 = (teta/n)*sv*(i1+phi*i21)-nu*v1
    dv2 = (teta/n)*sv*(i2+phi*i12)-nu*v2

    out = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2]

    return out

def rk(inp,t,dt,n):
    k1=field(inp,t,n)
    k2=field(soma_lista(inp,multiplica_lista((dt/2),k1)),t+dt/2,n)
    k3=field(soma_lista(inp,multiplica_lista((dt/2),k2)),t+dt/2,n)
    k4=field(soma_lista(inp,multiplica_lista(dt,k3)),t+dt,n)

    out = nega_zero(soma_lista(inp,multiplica_lista(dt/6,(soma_lista(soma_lista(k1,multiplica_lista(2,k2)),soma_lista(multiplica_lista(2,k3),k4))))))
    return out

def main():
    in0 = [None,None,None,None,None,None,None,None,None,None,None,None,None]

    in0[0]=700
    in0[1]=200
    in0[2]=100
    in0[3]=0
    in0[4]=0
    in0[5]=0
    in0[6]=0
    in0[7]=0
    in0[8]=0
    in0[9]=0
    in0[10]=9000
    in0[11]=500
    in0[12]=500


    t = 0

    dt = 1/365

    acum = 0

    for i in range(10):
        acum += in0[i]

    n0 = acum #numero de pessoas

    kk = 100*365

    """for k in range(kk):
        in0 = rk(in0,t,dt,n0)
        t = t+dt
        Para kk passos
"""

    in0 = rk(in0,t,dt,n0)
    print(in0)

main()

#print("digite os valores do input:")
#for i in range(13):
#print(f"valor {i+1}: ")
#inp.append(float(input()))
#print(f"lista do input: {inp}")
