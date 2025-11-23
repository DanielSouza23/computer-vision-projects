# Dolly Zoom Rendering & Stereo Vision Reconstruction

This project implements a 3D rendering–based dolly zoom effect, stereo image generation, and disparity map estimation using rank/census transforms and Hamming distance. 
It uses a point-cloud scene, custom camera poses, and Open3D for synthetic image creation.

Note: The dataset file `data/data.obj` (~267 MB) is not included in this repository due to its size.

---

## Features

### Dolly Zoom Rendering
- Loads foreground/background point clouds.
- Renders 75+ frames with changing focal length.
- Produces a smooth video (`dolly_zoom.mp4`).

### Stereo Image Generation
- Computes left/right images using a virtual stereo baseline.
- Converts rendered images to grayscale.

### Disparity Estimation
- Rank and Census transforms (5×5)
- Hamming distance cost computation
- Winner-Take-All matching
- Cost aggregation windows: 5×5, 9×9, 15×15

---

