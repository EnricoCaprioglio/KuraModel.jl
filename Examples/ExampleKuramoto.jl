# Example Kuramoto
using Kuramodel
greet()

Kura_step(
        [1,1],
        [1,1],
        [0 1; 1 0],
        0.02,
        θ=[2,0.5]
    )