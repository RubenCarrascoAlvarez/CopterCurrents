# DISCO
DISpersion relation based Currents Observations

                   --------------------------------------------------------
                   DISpersion relation based Currents Observations (DISCO)
                      	     1.00 Development Branch
                   ---------------------------------------------------------

This is the first stable release: DISCO 1.00, for detailed 
installation instructions, see the file INSTALL.



1. DISCO <br />
================

DISCO is a set of tools designed to extract water surface 
current fields (flow speeds) from UAV videos.

The Disco current extraction is based in the wave disperion relation, 
and the complete process is described in following work:

M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
doi: 10.1109/LGRS.2017.2749120
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8049342&isnumber=8082817

The DISCO signal processing is dividing in 3 steps:

	1) Record a nadir video of the water surface by a  
	quadcopter with a camera gimbal.

	2) Video rectification: Georeference the video in real world coordenates, 
	and slice the video in serveral squared regions. 

	3) Fit the most probable current, applying the wave dispersion relation 
	in a spectral energy-based maximization technique, in every squared region.



2. LICENSE <br />
================

* The DISCO application core, and other portions of the official DISCO
  distribution not explicitly licensed otherwise, are licensed under
  the GNU GENERAL PUBLIC LICENSE -- see the 'LICENSE' file in this
  directory for details.

* If you create a program which invokes (or provides) methods within
  (or for) the GPL DISCO application core, then the DISCO developers'
  position is that this is a 'mere aggregation' of the program 
  invoking the method and the program implementing the method as per
  section 2 of the GNU General Public License.




3. WEB RESOURCES (TO DO) <br />
================

DISCO's home page is at:

	https://www.xxxx.hzg/

Please be sure to visit this site for information, documentation,
tutorials, news, etc.

The latest version of DISCO can be found at:

	https://www.xxxx.hzg/downloads/



Have fun,

Ruben Carrasco <br />
Michael Stresser <br />
Jochen Horstmann <br />
