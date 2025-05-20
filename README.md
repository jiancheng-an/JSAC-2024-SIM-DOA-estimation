# JSAC-2024-SIM-DOA-estimation

This code provides the source files and generation codes of all the figure in this paper:
* J. An et al., "Two-Dimensional Direction-of-Arrival Estimation Using Stacked Intelligent Metasurfaces," in IEEE Journal on Selected Areas in Communications, vol. 42, no. 10, pp. 2786-2802, Oct. 2024.

# bib
@ARTICLE{JSAC_2024_An_Two,
  author={An, Jiancheng and Yuen, Chau and Guan, Yong Liang and Renzo, Marco Di and Debbah, Mérouane and Poor, H. Vincent and Hanzo, Lajos},
  journal={IEEE J. Sel. Areas Commun.}, 
  title={Two-Dimensional Direction-of-Arrival Estimation Using Stacked Intelligent Metasurfaces}, 
  year={2024},
  month={Oct.},
  volume={42},
  number={10},
  pages={2786-2802},
  doi={10.1109/JSAC.2024.3414613}}

# Figure List
* Fig. 1. Photographs of some existing SIMs. (a) A three-layer SIM; (b) A five-layer SIM.
* Fig. 2. The considered SIM-aided array system.
* Fig. 3. Top and front views of the considered SIM-based systems.
* Fig. 4. The proposed SIM-based DOA estimation protocol, where different colors refer to non-overlapping spatial frequency bins associated with different snapshots. This is achieved by adjusting the phase shifts of the input metasurface layer. The vertical and horizontal lines correspond to $\psi_x$ and $\psi_y$, respectively.
* Tab. I. Normalized loss function (in dB) based on a (2, 2) grid points
* Tab. II. Normalized loss function (in dB) based on a (4, 4) grid points
* Fig. 5 Convergence behavior of the proposed gradient descent algorithm based on a (2, 2) grid points
  * Fig. 5(a) Normalized loss function $\mathcal{L}$ versus the number of iterations
  * Fig. 5(b) Normalized loss function $\mathcal{L}$ versus the decay parameter $\zeta$