import matplotlib.pyplot as plt
import numpy as np

# Set seed for reproducibility
np.random.seed(993) # concordant
T1,T2=20,300

np.random.seed(99) # discordant
T1,T2=20,40

# Time parameters
T = 500
dt = 0.1
t_end = int(T / dt)
time = np.linspace(0, T, t_end)
T1,T2=20,40
# T1,T2=1,20

# Fix the index error by computing path_A up to t_end before trying to access index at T=4
# Simulate Brownian motion until T=2
t_bif1 = int(T1 / dt)
segment_ABC = np.cumsum( dt* np.random.randn(t_bif1))

# From T=2 to T=3, create two branches: A and B
t_bif2 = int(T2 / dt)
segment_AB = np.cumsum(dt * np.random.randn(t_bif2 - t_bif1)) + segment_ABC[-1]

# From T=3 to T=6, continue A and split B into B and C
segment_A = np.cumsum(dt * np.random.randn(t_end - t_bif2))+ segment_AB[-1] 
segment_B = np.cumsum(dt * np.random.randn(t_end - t_bif2))+ segment_AB[-1] 
segment_C = np.cumsum(dt * np.random.randn(t_end - t_bif1)) + segment_ABC[-1]

path_A = np.concatenate([segment_ABC,   segment_AB,  segment_A])
# path_B = np.concatenate([segment_ABC, segment_ABC[-1] + segment_AB, segment_AB[-1] + segment_B])
# path_C = np.concatenate([segment_ABC, segment_ABC[-1] + segment_C])
# Plotting
plt.figure(figsize=(10, 6))

plt.plot(time,path_A, label='Path A', color='blue') # A
plt.plot(time[time<=T1], segment_ABC,  color='black') # ABC

plt.plot(time[(T1<time) &( time<=T2)], segment_AB,  color='purple') # AB

plt.plot(time[time>T2], segment_B, label='Path B', color='red') # B

plt.plot(time[time>T1],  segment_C, label='Path C', color='green') # C

plt.axvline(x=T1, color='gray', linestyle='--', linewidth=1)
plt.axvline(x=T2, color='gray', linestyle='--', linewidth=1)
# plt.axvline(x=4, color='gray', linestyle='--', linewidth=1)
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)
plt.text(T1,-.01, 'T1', fontsize=12, color='black', transform=plt.gca().get_xaxis_transform(),
            ha='center', va='top')
plt.text(T2,-.01, 'T2', fontsize=12, color='black',transform=plt.gca().get_xaxis_transform(),
            ha='center', va='top')
plt.xlabel('Time')
plt.ylabel('Trait Value')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.xlim(0, T)
plt.savefig('figures/discordant.png')
plt.show()
 