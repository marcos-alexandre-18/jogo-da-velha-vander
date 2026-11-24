algoritmo "polinomios_com_coeficientes_reais"
// ------------------------------------------------------------------
// POLINOMIOS COM COEFICIENTES REAIS
//   p(x) = x^2 + 3x - 10
//   q(x) = x^3 - 3x^2 - 6x + 8
//   r(x) = 2x^3 - 10x^2 + 6x + 18
//   s(x) = x^4 - 10x^2 + 9
//
// a) raizes de p, q, r e s
// b) f = p * q  -> polinomio resultante e raizes
// c) g = p * r  -> polinomio resultante e raizes
// d) h = p * s  -> polinomio resultante e raizes
//
// As raizes sao calculadas pelo metodo de Durand-Kerner, que acha
// todas as raizes (reais e complexas) ao mesmo tempo.
// A[k, i] = coeficiente de x^i do polinomio k
//   k = 1 (p), 2 (q), 3 (r), 4 (s), 5 (f), 6 (g), 7 (h)
// ------------------------------------------------------------------
var
   A: vetor [1..7, 0..8] de real
   G: vetor [1..7] de inteiro
   NOME: vetor [1..7] de caractere
   C: vetor [0..8] de real
   RR, RI: vetor [1..8] de real
   i, k, op: inteiro
   x, v: real

// ---------- arredonda para 4 casas decimais ----------
funcao Arred(valor: real): real
var
   res: real
inicio
   se valor >= 0 entao
      res <- int(valor * 10000 + 0.5) / 10000
   senao
      res <- -int(-valor * 10000 + 0.5) / 10000
   fimse
   se abs(res) < 0.00005 entao
      res <- 0
   fimse
   retorne res
fimfuncao

// ---------- valor do polinomio k em um ponto (Horner) ----------
funcao Avalia(kk: inteiro; ponto: real): real
var
   m: inteiro
   soma: real
inicio
   soma <- 0
   para m de G[kk] ate 0 passo -1 faca
      soma <- soma * ponto + A[kk, m]
   fimpara
   retorne soma
fimfuncao

// ---------- le grau e coeficientes do polinomio k ----------
procedimento LerPolinomio(kk: inteiro)
var
   m: inteiro
inicio
   escreval("Polinomio ", NOME[kk], "(x):")
   repita
      escreva("  Grau (1 a 4): ")
      leia(G[kk])
   ate (G[kk] >= 1) e (G[kk] <= 4)
   para m de 0 ate G[kk] faca
      escreva("  Coeficiente de x^", m, ": ")
      leia(A[kk, m])
   fimpara
fimprocedimento

// ---------- multiplica polinomio u por v e guarda em w ----------
procedimento Multiplica(u, v, w: inteiro)
var
   m, n2: inteiro
inicio
   G[w] <- G[u] + G[v]
   para m de 0 ate 8 faca
      A[w, m] <- 0
   fimpara
   para m de 0 ate G[u] faca
      para n2 de 0 ate G[v] faca
         A[w, m + n2] <- A[w, m + n2] + A[u, m] * A[v, n2]
      fimpara
   fimpara
fimprocedimento

// ---------- escreve o polinomio k na tela ----------
procedimento MostraPolinomio(kk: inteiro)
var
   m: inteiro
   cf: real
   primeiro: logico
inicio
   escreva("  ", NOME[kk], "(x) = ")
   primeiro <- verdadeiro
   para m de G[kk] ate 0 passo -1 faca
      cf <- A[kk, m]
      se abs(cf) > 0.000001 entao
         // sinal
         se cf < 0 entao
            se primeiro entao
               escreva("-")
            senao
               escreva(" - ")
            fimse
         senao
            se nao primeiro entao
               escreva(" + ")
            fimse
         fimse
         // coeficiente (omite o 1 quando acompanha x)
         se (abs(abs(cf) - 1) > 0.000001) ou (m = 0) entao
            se abs(abs(cf) - int(abs(cf))) < 0.000001 entao
               escreva(int(abs(cf)))
            senao
               escreva(abs(cf):1:2)
            fimse
            se m > 0 entao
               escreva("*")
            fimse
         fimse
         // variavel
         se m = 1 entao
            escreva("x")
         fimse
         se m > 1 entao
            escreva("x^", m)
         fimse
         primeiro <- falso
      fimse
   fimpara
   escreval("")
fimprocedimento

// ---------- calcula e mostra as raizes do polinomio k ----------
procedimento Raizes(kk: inteiro)
var
   n2, m, j2, it, aux: inteiro
   raio, dr, di, pr, pim, nr, ni, t, den, ang, piv: real
inicio
   piv <- 3.14159265358979
   n2 <- G[kk]
   // descarta coeficientes lideres nulos
   enquanto (n2 > 0) e (abs(A[kk, n2]) < 0.000000001) faca
      n2 <- n2 - 1
   fimenquanto

   se n2 = 0 entao
      escreval("  Polinomio constante: nao tem raizes.")
   senao
      // normaliza (coeficiente lider = 1)
      para m de 0 ate n2 faca
         C[m] <- A[kk, m] / A[kk, n2]
      fimpara
      // chute inicial: pontos num circulo que contem todas as raizes
      raio <- 0
      para m de 0 ate n2 - 1 faca
         se abs(C[m]) > raio entao
            raio <- abs(C[m])
         fimse
      fimpara
      raio <- 1 + raio
      para m de 1 ate n2 faca
         ang <- 2 * piv * m / n2 + 0.4
         RR[m] <- raio * cos(ang)
         RI[m] <- raio * sen(ang)
      fimpara

      // iteracoes de Durand-Kerner
      para it de 1 ate 1000 faca
         para m de 1 ate n2 faca
            // p(z) por Horner, com z = RR[m] + i*RI[m]
            pr <- 1
            pim <- 0
            para j2 de n2 - 1 ate 0 passo -1 faca
               t <- pr * RR[m] - pim * RI[m] + C[j2]
               pim <- pr * RI[m] + pim * RR[m]
               pr <- t
            fimpara
            // denominador = produto de (z_m - z_j)
            dr <- 1
            di <- 0
            para j2 de 1 ate n2 faca
               se j2 <> m entao
                  nr <- RR[m] - RR[j2]
                  ni <- RI[m] - RI[j2]
                  t <- dr * nr - di * ni
                  di <- dr * ni + di * nr
                  dr <- t
               fimse
            fimpara
            den <- dr * dr + di * di
            se den > 0 entao
               nr <- (pr * dr + pim * di) / den
               ni <- (pim * dr - pr * di) / den
               RR[m] <- RR[m] - nr
               RI[m] <- RI[m] - ni
            fimse
         fimpara
      fimpara

      // ordena pela parte real (crescente)
      para m de 1 ate n2 - 1 faca
         para j2 de m + 1 ate n2 faca
            se RR[j2] < RR[m] entao
               t <- RR[m]
               RR[m] <- RR[j2]
               RR[j2] <- t
               t <- RI[m]
               RI[m] <- RI[j2]
               RI[j2] <- t
            fimse
         fimpara
      fimpara

      // mostra as raizes
      escreval("  Raizes de ", NOME[kk], "(x):")
      para m de 1 ate n2 faca
         se abs(RI[m]) < 0.0001 entao
            escreval("    x", m, " = ", Arred(RR[m]):1:4)
         senao
            se RI[m] >= 0 entao
               escreval("    x", m, " = ", Arred(RR[m]):1:4, " + ", Arred(RI[m]):1:4, "i")
            senao
               escreval("    x", m, " = ", Arred(RR[m]):1:4, " - ", Arred(-RI[m]):1:4, "i")
            fimse
         fimse
      fimpara
   fimse
fimprocedimento

inicio
   NOME[1] <- "p"
   NOME[2] <- "q"
   NOME[3] <- "r"
   NOME[4] <- "s"
   NOME[5] <- "f"
   NOME[6] <- "g"
   NOME[7] <- "h"

   // zera tudo
   para k de 1 ate 7 faca
      G[k] <- 0
      para i de 0 ate 8 faca
         A[k, i] <- 0
      fimpara
   fimpara

   escreval("=== POLINOMIOS COM COEFICIENTES REAIS ===")
   escreval("1 - Usar os polinomios da prova")
   escreval("2 - Digitar os coeficientes de p, q, r e s")
   escreva("Opcao: ")
   leia(op)
   escreval("")

   se op = 1 entao
      // p(x) = x^2 + 3x - 10
      G[1] <- 2
      A[1, 0] <- -10
      A[1, 1] <- 3
      A[1, 2] <- 1
      // q(x) = x^3 - 3x^2 - 6x + 8
      G[2] <- 3
      A[2, 0] <- 8
      A[2, 1] <- -6
      A[2, 2] <- -3
      A[2, 3] <- 1
      // r(x) = 2x^3 - 10x^2 + 6x + 18
      G[3] <- 3
      A[3, 0] <- 18
      A[3, 1] <- 6
      A[3, 2] <- -10
      A[3, 3] <- 2
      // s(x) = x^4 - 10x^2 + 9
      G[4] <- 4
      A[4, 0] <- 9
      A[4, 1] <- 0
      A[4, 2] <- -10
      A[4, 3] <- 0
      A[4, 4] <- 1
   senao
      para k de 1 ate 4 faca
         LerPolinomio(k)
         escreval("")
      fimpara
   fimse

   // produtos: f = p*q, g = p*r, h = p*s
   Multiplica(1, 2, 5)
   Multiplica(1, 3, 6)
   Multiplica(1, 4, 7)

   escreval("==================================================")
   escreval("a) RAIZES DE CADA POLINOMIO")
   escreval("==================================================")
   para k de 1 ate 4 faca
      MostraPolinomio(k)
      Raizes(k)
      escreval("")
   fimpara

   escreval("==================================================")
   escreval("b) f = p * q")
   escreval("==================================================")
   MostraPolinomio(5)
   Raizes(5)
   escreval("")

   escreval("==================================================")
   escreval("c) g = p * r")
   escreval("==================================================")
   MostraPolinomio(6)
   Raizes(6)
   escreval("")

   escreval("==================================================")
   escreval("d) h = p * s")
   escreval("==================================================")
   MostraPolinomio(7)
   Raizes(7)
   escreval("")

   // extra: valor dos polinomios para um x escolhido
   escreva("Calcular os polinomios para um valor de x? (1=sim, 0=nao): ")
   leia(op)
   enquanto op = 1 faca
      escreva("x = ")
      leia(x)
      para k de 1 ate 7 faca
         v <- Avalia(k, x)
         escreval("  ", NOME[k], "(", x, ") = ", v:1:4)
      fimpara
      escreva("Calcular outro x? (1=sim, 0=nao): ")
      leia(op)
   fimenquanto
fimalgoritmo
