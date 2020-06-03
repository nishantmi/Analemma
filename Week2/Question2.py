import astropy.units as u
from astropy.coordinates import SkyCoord

star1 = SkyCoord('00h08m',"+29d05m")
star2 = SkyCoord('23h04m',"+28d05m")
star3 = SkyCoord('23h05m',"+15d12m")
star4 = SkyCoord('00h13m',"+15d11m")

diag1 = star1.separation(star3)
diag2 = star2.separation(star4)

print(diag1)
print(diag2)
