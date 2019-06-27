function [pts] = rotatePoint3d(pts,yaw,pitch,roll)
% Rotate points according to drone yaw,pitch,roll
%
% https://en.wikipedia.org/wiki/Euler_angles
% https://de.wikipedia.org/wiki/Eulersche_Winkel
% Drehfolgen in der Fahrzeugtechnik
% Koordinatentransformation in der Flugtechnik, Drehfolge z, y′, x″ (Gier-Nick-Roll)
% blau: raumfestes Koordinatensystem (Ausgangslage)
% grün: Zwischenlage der y-Achse (Knotenachse N(y′))
% rot: körperfestes Koordinatensystem (Ziellage)
% Anmerkung: z-Achse zeigt in Praxis der Flugtechnik nach unten zur Erde (s. nächstes Bild)
% Gier-, Nick- und Rollwinkel(\Psi, \Theta, \Phi) als Lagewinkel eines Flugzeugs
% Gier- und Nickwinkel gegenüber erdfestem System,
% Rollwinkel gegenüber Flugzeug-Längsachse
% 
% Die in der Fahrzeugtechnik angewendeten und genormten (Luftfahrt: DIN 9300; Automobilbau: 
% DIN 70000; Schifffahrt) Drehfolgen gehören in die Gruppe der Tait-Bryan-Drehungen 
% (Drehung um die drei Koordinatenachsen). In den Normen sind die Verwendung der Formelzeichen
% \Psi, \Theta und \Phi und die Namen Gier-, Nick- und Roll-Winkel (engl. yaw, pitch and roll 
% angle) für die drei Eulerwinkel vorgeschrieben.[Anmerkungen 3] Durch die drei Drehungen wird
% bdas erdfeste System (engl. world frame) xyz in das körperfeste Koordinatensystem 
% (engl. body frame) XYZ oder umgekehrt gedreht.
% Gier-Nick-Roll: z, y′, x″-Konvention
% 
% Die Drehfolge ist:
% 
%     Mit dem im erdfesten System gemessenen Gierwinkel \Psi (auch Steuerkurs oder Azimut genannt) wird um die z-Achse gedreht. Die y-Achse wird zur Knotenachse N(y').
%     Hauptwertebereich: - \pi < \Psi \le \pi
%     Mit dem gegen die Erdoberfläche (x-y-Ebene) gemessenen Nickwinkel \Theta wird um die Knotenachse N(y') gedreht. Es entsteht die fahrzeugfeste X-Achse.
%     Hauptwertebereich: - \frac{\pi}{2} \le \Theta \le \frac{\pi}{2}
%     Der Rollwinkel \Phi (auch Wankwinkel genannt) beschreibt die Drehung um die fahrzeugfeste X-Achse. Es entstehen die fahrzeugfesten Achsen Y und Z.
%     Hauptwertebereich: - \pi < \Phi \le \pi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function coded by Michael Stresser for CopterCurrents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of CopterCurrents
%
% CopterCurrents has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% CopterCurrents is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CopterCurrents is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CopterCurrents.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rr = [1 0 0; 0 cosd(roll) sind(roll); 0 -sind(roll) cosd(roll)];
Rp = [cosd(pitch) 0 -sind(pitch); 0 1 0; sind(pitch) 0 cosd(pitch)];
Ry = [cosd(yaw) sind(yaw) 0; -sind(yaw) cosd(yaw) 0; 0 0 1];

Rt = Rr*Rp*Ry;

pts = Rt * pts';
pts=pts';

end
