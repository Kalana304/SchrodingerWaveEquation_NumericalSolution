# Numerical Solution for 1-D Schrödinger’s Wave Equation
Quantum Physics, unlike the Classical Physics, describes the behavior of the matter and light (photons) on the atomic and subatomic level. It is heavily based on the mathematical models which attempt to describe and account for the properties of atoms and their constituents. Schrödinger Equation, also known as Schrödinger’s Wave Equation, is a partial differential equation that describes the dynamics of a quantum mechanics system using the wave function.

In this project, a numerical solution has been proposed to solve the 1-D Schrödinger’s Wave Equation for a user defined Energy Function using python. Further the idea has been extended to 2-D Wave Functions as well. 

The Time-Dependent Schrödinger Equation is as follows: −ħ22𝑚∇2𝛹(𝑟⃗,𝑡)+𝑉(𝑟⃗)𝛹(𝑟⃗,𝑡)=𝑖ħ𝜕𝛹(𝑟⃗,𝑡)𝜕𝑡
Where,
ħ−𝑟𝑒𝑑𝑢𝑐𝑒𝑑 𝑃𝑙𝑎𝑛𝑐𝑘 𝑐𝑜𝑛𝑠𝑡𝑎𝑛𝑡 (1.0546 × 10−34𝐽𝑠)
𝑚−𝐸𝑓𝑓𝑒𝑐𝑡𝑖𝑣𝑒 𝑚𝑎𝑠𝑠 𝑜𝑓 𝑡ℎ𝑒 𝑞𝑢𝑎𝑛𝑡𝑢𝑚 𝑝𝑎𝑟𝑡𝑖𝑐𝑙𝑒
𝛹(𝑟⃗,𝑡)−𝑇𝑖𝑚𝑒 𝐷𝑒𝑝𝑒𝑛𝑑𝑒𝑛𝑡 𝑤𝑎𝑣𝑒 𝑓𝑢𝑛𝑐𝑡𝑖𝑜𝑛
𝑟⃗−𝑃𝑜𝑠𝑖𝑡𝑖𝑜𝑛 𝑉𝑒𝑐𝑡𝑜𝑟 𝑜𝑓 𝑡ℎ𝑒 𝑝𝑎𝑟𝑡𝑖𝑐𝑙𝑒
𝑉(𝑟⃗)−𝑃𝑜𝑡𝑒𝑛𝑠𝑡𝑖𝑎𝑙 𝐸𝑛𝑒𝑟𝑔𝑦 𝐹𝑢𝑛𝑐𝑡𝑖𝑜𝑛
Since we are more interested in stationary quantum states rather than how they evolve with time, we usually consider the Time-Independent version of the Schrödinger’s wave equation as below: −ħ22𝑚∇2𝛹(𝑟⃗)+𝑉(𝑟⃗)𝛹(𝑟⃗)=𝐸𝛹(𝑟⃗)
Where,
𝛹(𝑟⃗)−𝑇𝑖𝑚𝑒 𝐼𝑛𝑑𝑒𝑝𝑒𝑛𝑑𝑒𝑛𝑡 𝑤𝑎𝑣𝑒 𝑓𝑢𝑛𝑐𝑡𝑖𝑜𝑛
𝐸−𝑇𝑜𝑡𝑎𝑙 𝑒𝑛𝑒𝑟𝑔𝑦 𝑜𝑓 𝑡ℎ𝑒 𝑝𝑎𝑟𝑡𝑖𝑐𝑙𝑒
According to the Heisenberg’s Uncertainty Principle, we cannot predict the exact position and momentum of a particle at the same time in quantum world. We can only predict a probability for the particle being at a specific location 𝑟⃗. By solving the Schrödinger equation, we can get the Probability Density Function (PDF) for a particle being at 𝑟⃗ as follows.
𝑝(𝑟⃗)= |𝛹(𝑟⃗)|2= 𝛹(𝑟⃗)×𝛹(𝑟⃗)∗
In this report, we explore a method to numerically solve 1-Dimensional wave equation for a user specified 1-Dimensional Potential Energy Function 𝑉(𝑥) using Python. The solution will calculate the user defined

