import netCDF4
import numpy as np
from numpy.lib.recfunctions import append_fields


def main():
    nc_fname = 'NISKINE-MR-Near-Surface-Currents.nc'
#    nc_fname = '/path/to/NISKINE-MR-Near-Surface-Currents.nc'
    grid_res = 500.
    max_dt = 11.1 / (24. * 60.)    

    with netCDF4.Dataset(nc_fname) as src:
        ux = src['ux'][:]
        uy = src['uy'][:]
#        depth = src['depth'][:]
        start_time = src['startTime'][:]
        end_time = src['endTime'][:]
        lat = src['latitude'][:]
        lon = src['longitude'][:]
        mr = np.core.records.fromarrays(
            [0.5 * (start_time + end_time), lon, lat, ux, uy],
            names='y_day, lon, lat, ux, uy'
        )

    mr = compute_flow_characteristics(mr, grid_res, max_dt)
    print "Done in main"
    
def compute_flow_characteristics(mr, grid_res, max_dt):

    print "in compute..."

    vorticity = np.full(mr.size, np.nan)
    divergence = np.full(mr.size, np.nan)

    print mr.size

    for index, mr_select in enumerate(mr):

        print (1.*index/mr.size)*100,'% done'

        foo = np.where((mr.y_day > (mr_select.y_day - max_dt)) & (mr.y_day < (mr_select.y_day + max_dt)))[0]
        mr_ = mr[foo]
        foo = flat_earth(mr_select.lon, mr_select.lat, lat1=mr_.lat, lon1=mr_.lon)
        mr_ = append_fields(mr_, ('x', 'y'), (foo.dist_e, foo.dist_n), asrecarray=True, usemask=False)
        dx = np.round(mr_.x / grid_res).astype(int)
        dy = np.round(mr_.y / grid_res).astype(int)
        neighbor_id = np.where((np.abs(dx) <= 1) & (np.abs(dy) <= 1) & ((np.abs(dx) + np.abs(dy)) != 0))[0]
        if neighbor_id.size == 0:
            continue
        mr_ = mr_[neighbor_id]
        dx = dx[neighbor_id]
        dy = dy[neighbor_id]
        left_id = np.where((dx < 0) & (dy == 0))[0]
        right_id = np.where((dx > 0) & (dy == 0))[0]
        bottom_id = np.where((dx == 0) & (dy < 0))[0]
        top_id = np.where((dx == 0) & (dy > 0))[0]
        if ((left_id.size == 0) | (right_id.size == 0) | (bottom_id.size == 0) | (top_id.size == 0)):
            continue
        if left_id.size > 1:
            left_id = left_id[np.argmin(np.abs(mr_.y_day[left_id] - mr_select.y_day))]
        else:
            left_id = left_id[0]
        if right_id.size > 1:
            right_id = right_id[np.argmin(np.abs(mr_.y_day[right_id] - mr_select.y_day))]
        else:
            right_id = right_id[0]
        if bottom_id.size > 1:
            bottom_id = bottom_id[np.argmin(np.abs(mr_.y_day[bottom_id] - mr_select.y_day))]
        else:
            bottom_id = bottom_id[0]
        if top_id.size > 1:
            top_id = top_id[np.argmin(np.abs(mr_.y_day[top_id] - mr_select.y_day))]
        else:
            top_id = top_id[0]
        dx = 0.5 * (mr_.x[right_id] - mr_.x[left_id])
        dy = 0.5 * (mr_.y[top_id] - mr_.y[bottom_id])
        divergence[index] = compute_divergence(
            mr_.ux[left_id], mr_.ux[right_id], dx,
            mr_.uy[bottom_id], mr_.uy[top_id], dy,
            lat=mr_select.lat
        )
        vorticity[index] = compute_vorticity(
            mr_.ux[bottom_id], mr_.ux[top_id], dx,
            mr_.uy[left_id], mr_.uy[right_id], dy,
            lat=mr_select.lat
        )
    mr = append_fields(mr, ('divergence', 'vorticity'), (divergence, vorticity), asrecarray=True, usemask=False)

    print "returning"

    return mr

def flat_earth(lon0, lat0, bearing=None, distance=None,
               distance_east=None, distance_north=None,
               lat1=None, lon1=None):

    """Flat earth approximation: From reference position (lon0, lat0) and
    a second positon in the vicinity thereof to distance and bearing
    (and vice versa).

    Approximation fails in the vicinity of either
    pole and at large distances.

    Fractional errors are of order (distance/R)**2.

    To compare this solution with an excact one from the geographiclib
    module, execute the following lines of Python code:
    import math
    from geographiclib.geodesic import Geodesic
    geod = Geodesic.WGS84
    g = geod.Inverse(lat0, lon0, lat1, lon1)
    print "The distance is {:.3f} m.".format(g['s12'])
    print "The azimuth is {:.3f} m.".format((450.-g['azi1']) % 360.)
    g = geod.Direct(lat0, lon0, bearing, distance)
    print "The position is ({:.8f}, {:.8f}).".format(g['lat2'],g['lon2'])
    
    Source for the following equations:
    http://edwilliams.org/avform.htm (linked to by
    https://www.nhc.noaa.gov)

    Args:

    lon0 (float) - Reference longitude in decimal degrees

    lat0 (float) - Reference latitude in decimal degrees

    bearing (float, optional) - Bearing with reference position as
                                origin in deg clockwise from north;
                                default is None

    distance (float, optional) - Distance from reference position in
                                 meters; default is None

    distance_east(float, optional) - Distance along east-west axis
                                     from reference position in
                                     meters; default is None

    distance_north(float, optional) - Distance along north-south axis
                                      from reference position in
                                      meters; default is None

    lon1 (float, optional) - Longitude of second position in vicinity
                             to the reference position in decimal
                             degrees; default is None

    lat1 (float, optional) - Latitude of second position in vicinity
                             to the reference position in decimal
                             degrees; default is None

    """

    f = 1./298.257223563 # Flattening (WGS 84)
    a = 6378137. # [m] Equatorial radius of Earth (WGS 84)
    e_squared = f*(2-f)
    lat0_rad = np.radians(lat0)
    r1 = (a*(1.-e_squared)) / (1.-e_squared*np.sin(lat0_rad)**2)**(3./2.) # Meridional radius of curvature
    r2 = a / np.sqrt(1.-e_squared*np.sin(lat0_rad)**2) # Radius of curvarture in the prime vertical

    class ReturnValue(object):
        pass
        
    def inverse(lon1, lat1):
        d_lat = lat1-lat0
        d_lon = lon1-lon0
        dist_n = r1*np.radians(d_lat)
        dist_e = r2*np.cos(lat0_rad)*np.radians(d_lon)
        dist = np.sqrt(dist_n**2 + dist_e**2)
        bearing = np.degrees(np.arctan2(dist_n, dist_e))
        result = ReturnValue()
        result.dist_n = dist_n
        result.dist_e = dist_e
        result.dist = dist
        result.bearing = bearing
        return result
        
    def direct(dist_e, dist_n):
        d_lat_rad = dist_n/r1
        d_lon_rad = dist_e/(r2*np.cos(lat0_rad))
        result = ReturnValue()
        result.lat = lat0+np.degrees(d_lat_rad)
        result.lon = lon0+np.degrees(d_lon_rad)
        return result

    if lon1 is not None and lat1 is not None:
        return inverse(lon1, lat1)

    if distance_east is not None and distance_north is not None:
        return direct(distance_east, distance_north)
       
    if bearing is not None and distance is not None:
        distance_east = distance*np.cos(np.radians(bearing))
        distance_north = distance*np.sin(np.radians(bearing))
        return direct(distance_east, distance_north)

def compute_divergence(_ux, ux_, dx, _uy, uy_, dy, lat=None):
    """Compute normalized divergence (positive -> divergence (upwelling);
    negative -> convergence (downwelling); [1/s]) at one location
    within a 2D grid.

    Args:

    _ux (float) - The x component of the current located one grid cell to the left [m/s]

    ux_ (float) - The x component of the current located one grid cell to the right [m/s]

    dx (float) - The grid resolution along the x axis [m]

    _uy (float) - The y component of the current located one grid cell below [m/s]

    uy_ (float) - The y component of the current located one grid cell above [m/s]

    dy (float) - The grid resolution along the y axis [m]

    """
    divergence = ((ux_ - _ux) / (2. * dx) +
                  (uy_ - _uy) / (2. * dy))
    if lat is not None:
        divergence /= inertial_frequency(lat)
    return divergence

def compute_vorticity(_ux, ux_, dx, _uy, uy_, dy, lat=None):
    """Compute normalized vorticity (positive -> counter-clockwise;
    negative -> clockwise; [1/s]) at one location within a 2D grid.

    Args:

    _ux (float) - The x component of the current located one grid cell below [m/s]

    ux_ (float) - The x component of the current located one grid cell above [m/s]

    dx (float) - The grid resolution along the x axis [m]

    _uy (float) - The y component of the current located one grid cell to the left [m/s]

    uy_ (float) - The y component of the current located one grid cell to the right [m/s]

    dy (float) - The grid resolution along the y axis [m]

    """
    vorticity = ((uy_ - _uy) / (2. * dx) -
                 (ux_ - _ux) / (2. * dy))
    if lat is not None:
        vorticity /= inertial_frequency(lat)
    return vorticity

def inertial_frequency(lat, return_period=False):
    omega = 7.2921 * 10. ** (-5.) # Earth rotation rate [rad/s]
    inertial_frequency = 2. * omega * np.sin(np.radians(lat)) # [rad/s]
    if return_period:
        inertial_period = (2. * np.pi / inertial_frequency) / (60. * 60.) # [h]
        return inertial_frequency, inertial_period
    return inertial_frequency

if __name__ == '__main__':
    sys.exit(main())
