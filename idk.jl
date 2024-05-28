# not sure what's happening below...

J = 1000
nxt = rand(Multinomial(J,fill(1/J,J)))
nxtₙ= rand(Multinomial(J,nxt/J))
nxt3= rand(Multinomial(J,nxtₙ/J))

var(nxt)

var(nxtₙ)

var(nxt3)


rand(Binomial(1000000,0.5))

thing1 = randn((2,3))

thing2 = rand(Bernoulli(0.5),3)

thing1*thing2

dot(thing1[1,:],thing2)
