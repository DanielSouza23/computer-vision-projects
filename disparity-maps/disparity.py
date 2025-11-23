import numpy as np
import imageio.v2 as iio
import matplotlib.pyplot as plt

L = iio.imread("teddyL.pgm").astype(np.uint8)
R = iio.imread("teddyR.pgm").astype(np.uint8)
GT = iio.imread("disp2.pgm").astype(np.uint8)  

# Part 1: Rank transform 5x5
def rank5(img):
    r = 2 # 5x5, radius = floor(win/2)
    H, W = img.shape
    out = np.zeros((H, W), np.uint8) # Border    
    for y in range(r, H-r): # only where a 5×5 - this means that image will have a black border of 2 pixels (not applicable for 5x5 as it wont have enough neighbors)
        for x in range(r, W-r):
            center = img[y, x]
            win = img[y-r:y+r+1, x-r:x+r+1] # 5×5 window around pixel
            out[y, x] = np.sum(win < center)   # count how many neighbors are strictly less than center
    return out

L_rank = rank5(L)
R_rank = rank5(R)

# Part 2: Census
def census5(img):
    r = 2 # 5x5
    H, W = img.shape
    out = np.zeros((H, W), np.uint32)
    for y in range(r, H - r): # only where 5x5 fits
        for x in range(r, W - r):           
            code = 0
            for dy in range(-r, r + 1):    # top-left to bottom-right
                for dx in range(-r, r + 1):
                    if dy == 0 and dx == 0:
                        continue # skip center
                    code <<= 1 # make room for next bit (new LSB, shift order to the left)
                    if img[y + dy, x + dx] < img[y, x]:
                        code |= 1 # x-20, example (slide 55)
                    else:
                        code |= 0 # x+35, example 
            out[y, x] = code # output bit string, 1s and 0s, top left to bottom right            
    return out

L_cen = census5(L_rank)
R_cen = census5(R_rank)

# Part 3: Hamming 5x5, 9x9, 15x15
def hamming(L, R):
    return int((int(L) ^ int(R)).bit_count()) # Compare L and R, see where bit strings differ, and count how many differences

def cost(L_cen, R_cen, d): # Cost  for each disparity
    H, W = L_cen.shape
    cost = np.full((H, W), 24, np.uint8)  # initialize with max possible Hamming cost (24) for unaccounted pixels
    x_start = d  # d = xl - xr
    for y in range(H): # row search
        for x in range(x_start, W): # compare pixels along epipolar line (shifted by d) 
            cost[y, x] = hamming(L_cen[y, x], R_cen[y, x - d])
    return cost

def sums(img, win): # Aggregation (Slide 32)
    r = win // 2 # windows = 5, 9, 15
    H, W = img.shape
    out = np.zeros_like(img, dtype=np.int32) # array of 0s same size as img
    for y in range(H): 
        for x in range(W):
            # Create box around pixel
            t_ = (y-r) # tentative
            if t_ < 0:
                top = 0
            else:
                top = t_
            b_ = (y+r+1) # tentative
            if b_ > H:
                bottom = H
            else:
                bottom = b_
            r_ = (x+r+1) # tentative
            if r_ > W:
                right = W
            else:
                right = r_
            l_ = (x-r) # tentative
            if l_ < 0:
                left = 0
            else:
                left = l_
            out[y, x] = img[top:bottom, left:right].sum() # Sum costs within box
    return out

# Part 4: Disparity Maps

def WTA(L_cen, R_cen, max_disp, win):
    H, W = L_cen.shape
    best = np.full((H, W), np.iinfo(np.int32).max, dtype=np.int32)   # store smallest aggregated cost, initial high value to easily be replaced
    disp = np.zeros((H, W), np.uint8) # result disparity for single left image pixel

    for d in range(max_disp + 1): # 64 disparities, search along right-image epipolar line
        cost_d = cost(L_cen, R_cen, d) # cost
        agg_d  = sums(cost_d, win) # aggregated cost 
        
        for y in range(H):  
            for x in range(W):
                # Compare calculated cost with current best
                if agg_d[y,x] < best[y,x]:
                    best[y,x] = agg_d[y,x]
                    disp[y,x] = d

    return disp

disp5  = WTA(L_cen, R_cen, 63, 5)
disp9  = WTA(L_cen, R_cen, 63, 9)
disp15 = WTA(L_cen, R_cen, 63, 15)

# Plot disparity maps
fig, ax = plt.subplots(1,3, figsize=(14,4)) 
for i,(D,w) in enumerate(((disp5,5),(disp9,9),(disp15,15))):
    ax[i].imshow(D, cmap="gray", vmin=0, vmax=63)
    ax[i].set_title(f"Disparity ({w}x{w})"); ax[i].axis("off")
plt.tight_layout()
plt.show()
fig.savefig("disp.png")


# Part 5: GT Comparison
GT4 = np.rint(GT.astype(np.float32)/4).astype(np.uint8) # Normalize GT disparity to match scale (0-63) of computed maps 

def bad_rate(dist, gt):
    diff = np.abs(dist.astype(np.int16) - gt.astype(np.int16)) # Error
    bad = (diff > 1) # If error > tolerance, bad pixel (boolean)
    return 100*bad.mean() # (since it is a boolean, 100*mean = percentage)

for D, w in ((disp5, 5), (disp9, 9), (disp15, 15)):
    error = bad_rate(D, GT4)
    print(f"Bad-pixel rate ({w}x{w}): {error:.2f}%")

