import os
import numpy
import gdspy
import numpy as np


print('Using gdspy module version ' + gdspy.__version__)


## ------------------------------------------------------------------ ##
##	POLYGONS														  ##
## ------------------------------------------------------------------ ##


## First we need a cell to add the polygons to.
poly_cell = gdspy.Cell('POLYGONS')



def trenches(x0, y0, dx, dy, n, layer, stagger):
	xs = numpy.arange(x0, 2*n*dx + x0, 2*dx)
	for x in xs:
		if stagger and int((x-x0)/(2.*dx))%2==0:
			poly_cell.add(gdspy.Rectangle((x, y0+500.), (x+dx, y0 + dy), layer))
		else:
			poly_cell.add(gdspy.Rectangle((x, y0), (x+dx, y0 + dy), layer))
        poly_cell.add(gdspy.Rectangle((x0, y0), (x0 + (2*n - 1)*dx, y0 + dx), layer))
	
def trench_pattern(x0, y0, sw, sh, bw, bh, spread_length, n, layer):
	center1 = x0+(2*n - 1)*sw/2.
	center2 = x0+(2*n - 1)*bw/2.
	lx0 = x0-center2 + center1
	ly0 = y0-bh-spread_length
	trenches(x0, y0, sw, sh, n, layer, False)
	#trenches(lx0, ly0, bw, bh, n, layer, False)
	lxs = numpy.arange(lx0, 2*n*bw + lx0, 2*bw)
	xs = numpy.arange(x0, 2*n*sw + x0, 2*sw)
	#for i in range(len(lxs)):
		#points = [(lxs[i], ly0+bh), (xs[i], y0), (xs[i]+sw, y0), (lxs[i]+bw, ly0+bh)]
		#poly = gdspy.Polygon(points, layer)
		#poly_cell.add(poly)


def outline(x0, y0, bw, bh, sw, sh, top_layer, bottom_layer, shield = 1.):
	gap = (bw-sw)/2.
        tw = 50.
        y0 += shield
	#poly_cell.add(gdspy.Rectangle((x0, y0), (x0-tw, y0 + bh), top_layer))
	poly_cell.add(gdspy.Rectangle((x0-tw, y0+bh), (x0+gap, y0+bh+tw), top_layer))
	poly_cell.add(gdspy.Rectangle((x0+(bw - gap), y0+bh), (x0+bw+tw, y0+bh+tw), top_layer))
	#poly_cell.add(gdspy.Rectangle((x0+bw+tw, y0), (x0+bw, y0 + bh), top_layer))
	#poly_cell.add(gdspy.Rectangle((x0-tw, y0-tw), (x0+bw+tw, y0), top_layer))
	poly_cell.add(gdspy.Rectangle((x0+gap-tw, y0+bh+tw), (x0+gap, y0+bh+sh), top_layer))
	poly_cell.add(gdspy.Rectangle((x0+gap-tw, y0+bh+sh), (x0+bw-gap+tw, y0 + bh+sh+tw), top_layer))
	poly_cell.add(gdspy.Rectangle((x0+bw-gap, y0+bh+tw), (x0+bw-gap+tw, y0+bh+sh), top_layer))	

	poly_cell.add(gdspy.Rectangle((x0, y0+bh+1250.), (x0+bw, y0 + bh), bottom_layer))
        #poly_cell.add(gdspy.Rectangle((x0, y0), (x0-500, y0 + bh+1250), bottom_layer))
        #poly_cell.add(gdspy.Rectangle((x0+bw, y0), (x0+bw+500, y0 + bh+1250), bottom_layer))
        #poly_cell.add(gdspy.Rectangle((x0-500, y0+500), (x0+bw+500, y0-500), bottom_layer))


def cantilever(x0, y0, bw, bh, sw, sh, trench_width, bottom_layer, shield, n):
	ntrenches = int(sw/(2.*trench_width)) - 1
	bgap = (bw-sw)/2.
	lgap = (sw - trench_width*(2*ntrenches - 1))/2.    
	outline(x0, y0-600., bw, bh+600., sw, sh,  bottom_layer + 2, bottom_layer, shield = shield)
	wire_pos = bh-sh-50
	trench_pattern(x0+lgap+bgap, y0+bh+sh - 100., trench_width, 100., 75, wire_pos, 500, ntrenches, bottom_layer+2)
        #Connection between trenches and plating bar
        poly_cell.add(gdspy.Rectangle((x0+bw/2.-trench_width/2., y0-600.), (x0+bw/2.+trench_width/2., y0+bh+sh-100.), bottom_layer + 2))
        poly_cell.add(gdspy.Rectangle((x0+bw/2.-trench_width*8./2., y0-600.), (x0+bw/2.+trench_width*8./2., y0+bh-100.), bottom_layer + 2))
        #For layer where we want to keep the gold for electroforming.
        trench_pattern(x0+lgap+bgap, y0+bh+sh - 100., trench_width, 100., 75, wire_pos, 500, ntrenches, bottom_layer+1)
        poly_cell.add(gdspy.Rectangle((x0+bw/2.-trench_width/2., y0-600.), (x0+bw/2.+trench_width/2., y0+bh+sh-100.), bottom_layer + 1))
        poly_cell.add(gdspy.Rectangle((x0+bw/2.-trench_width*8./2., y0-600.), (x0+bw/2.+trench_width*8./2., y0+bh-100.), bottom_layer + 1))
        poly_cell.add(gdspy.Text(str(shield), 300, (x0+800, y0+400), layer = bottom_layer + 2))
        poly_cell.add(gdspy.Text(str(n), 300, (x0+2000, y0+400), layer = bottom_layer + 2))
        
        #poly_cell.add(gdspy.Rectangle((x0-500, y0+500), (x0+bw+500, y0-500), bottom_layer))
        
#	poly_cell.add(gdspy.Rectangle((x0, y0-450), (x0+1500, y0+bh), bottom_layer+2))
#	poly_cell.add(gdspy.Rectangle((x0+bw, y0-450), (x0+bw-1500, y0+bh), bottom_layer+2))


def shield_thickness(i, increments = np.array([0, 2, 4, 6, 8, 10]), shields = [0., 2., 4., 6., 8., 10.]):
        
        return shields[np.argmin(np.abs(i**2 - increments**2))]
        
        #for ind in increments:
        #        if ind>i:
        #                n+= 1.
        #        if n>1.:
        #                print 'yay'
        #                return shields[increments.index(ind)]
                        
        
        

def cantilever_row(x_spacing, nx, x_start, y, ny, y_spacing, y0, r):
	xs = numpy.arange(x_start, nx*x_spacing+x_start, x_spacing)
        shields = numpy.array([1, 2, 4, 6, 8, 10])
        tshields = shields
	for i in range(len(xs)):
                x = xs[i]
                 
                s = numpy.random.choice(tshields)
                tshields = tshields[tshields != s]
                
                if not tshields.any():
                        tshields = shields
		cantilever(x-575., y, 3000, 5000, 500, 500, 25, 1, s, i + r*nx)
                
                #cross(x, y-500, 100, 500, 3)
                #cross(x, y+5000, 100, 500, 3) 
                #cross(x+5000, y-500, 100, 500, 3 )
                #cross(x+5000, y+5000, 100, 500, 3 )
        
	poly_cell.add(gdspy.Rectangle((x_start-1000, y-600), (x_start+(nx-0.25)*x_spacing, y-400), 3))
        #Copy for mask layer to keep gold
        poly_cell.add(gdspy.Rectangle((x_start-1000, y-600), (x_start+(nx-0.25)*x_spacing, y-400), 2))
        
        


def cantilever_grid(x0, y0, nx, ny, x_spacing, y_spacing):
	ys = numpy.arange(y0, ny*y_spacing + y0, y_spacing)
        xs = numpy.arange(x0, nx*x_spacing + x0, x_spacing)
        wire_width = 200.
	for i in range(len(ys)):
		 cantilever_row(x_spacing, nx, x0, ys[i], ny, y_spacing, y0, i)

        for x in xs:
                poly_cell.add(gdspy.Rectangle((x-1050, y0-600), (x-850, (ny)*y_spacing + y0 - 400), 3))
                #copy to protect gold
                poly_cell.add(gdspy.Rectangle((x-1050, y0-600), (x-850, (ny)*y_spacing + y0 - 400), 2))
        #Thin wire connectors                 
	#poly_cell.add(gdspy.Rectangle((x0-1050, y0-600), (x0-850, (ny)*y_spacing + y0-400), 3))
        poly_cell.add(gdspy.Rectangle((x0-1050 + nx*x_spacing, y0-600), (x0-850 + nx*x_spacing, (ny)*y_spacing + y0-400), 3))
        poly_cell.add(gdspy.Rectangle((x0-1050, y0-600+ny*y_spacing), (x0+(nx-0.25)*x_spacing, y0-400+ny*y_spacing), 3))
        #copies to protect gold seed layer
        poly_cell.add(gdspy.Rectangle((x0-1050 + nx*x_spacing, y0-600), (x0-850 + nx*x_spacing, (ny)*y_spacing + y0-400), 2))
        poly_cell.add(gdspy.Rectangle((x0-1050, y0-600+ny*y_spacing), (x0+(nx-0.25)*x_spacing, y0-400+ny*y_spacing), 2))

        #Fat contact pads
        lx = 10000.
        ly = 10000.
        osetx = 1050.
        osety = 500.
        

	poly_cell.add(gdspy.Rectangle((x0-lx - osetx, y0+ny*y_spacing/2.-ly/2.), (x0 - osetx, y0+ny*y_spacing/2.+ly/2.), 3))

        poly_cell.add(gdspy.Rectangle((x0 + lx + nx*x_spacing - osetx + wire_width, y0+ny*y_spacing/2.-ly/2.), (x0 + nx*x_spacing - osetx + wire_width, y0+ny*y_spacing/2.+ly/2.), 3))

        poly_cell.add(gdspy.Rectangle((x0 + nx*x_spacing/2. - lx/2. - osetx - wire_width/2. + 12.5, y0 + ny*y_spacing - osety + wire_width/2.), (x0 + nx*x_spacing/2. + lx/2. - osetx - wire_width/2. + 12.5, y0 + ly + ny*y_spacing - osety + wire_width/2.), 3))

        poly_cell.add(gdspy.Rectangle((x0 + nx*x_spacing/2. - lx/2. - osetx - wire_width/2. + 12.5, y0 - osety - wire_width/2.), (x0 + nx*x_spacing/2. + lx/2. - osetx - wire_width/2. + 12.5, y0 - ly - osety - wire_width/2.), 3))

        #copies to protect gold seed layer
        poly_cell.add(gdspy.Rectangle((x0-lx - osetx, y0+ny*y_spacing/2.-ly/2.), (x0 - osetx, y0+ny*y_spacing/2.+ly/2.), 2))

        poly_cell.add(gdspy.Rectangle((x0 + lx + nx*x_spacing - osetx + wire_width, y0+ny*y_spacing/2.-ly/2.), (x0 + nx*x_spacing - osetx + wire_width, y0+ny*y_spacing/2.+ly/2.), 2))

        poly_cell.add(gdspy.Rectangle((x0 + nx*x_spacing/2. - lx/2. - osetx - wire_width/2. + 12.5, y0 + ny*y_spacing - osety + wire_width/2.), (x0 + nx*x_spacing/2. + lx/2. - osetx - wire_width/2. + 12.5, y0 + ly + ny*y_spacing - osety + wire_width/2.), 2))

        poly_cell.add(gdspy.Rectangle((x0 + nx*x_spacing/2. - lx/2. - osetx - wire_width/2. + 12.5, y0 - osety - wire_width/2.), (x0 + nx*x_spacing/2. + lx/2. - osetx - wire_width/2. + 12.5, y0 - ly - osety - wire_width/2.), 2))

        
        
        
        
        


def allignment_mark_box(x0, y0, layer, w = 10000, t = 1000):

        poly_cell.add(gdspy.Rectangle((x0-w/2., y0-w/2.), (x0-w/2.-t, y0+w/2), layer))  
        poly_cell.add(gdspy.Rectangle((x0-w/2., y0+w/2.), (x0+w/2., y0+w/2+t), layer))  
        poly_cell.add(gdspy.Rectangle((x0+w/2., y0-w/2.), (x0+w/2.+t, y0+w/2), layer))  
        poly_cell.add(gdspy.Rectangle((x0-w/2., y0-w/2.), (x0+w/2., y0-w/2-t), layer))  
        
def box_grid(x0, y0, w, s, n, layer):
        tw = (n*w + (n-1.)*s)
        xl = x0 - tw/2.
        yb = y0 - tw/2.
        xs = np.arange(xl, xl+tw, w+s)
        ys = np.arange(yb, yb+tw, w+s)
        for x in xs:
                for y in ys:
                        poly_cell.add(gdspy.Rectangle((x, y), (x+w, y+w), layer))  
        
def cross(x0, y0, t, l, layer):
        poly_cell.add(gdspy.Rectangle((x0-l/2., y0-t/2.), (x0+l/2., y0+t/2.), layer))  
        poly_cell.add(gdspy.Rectangle((x0-t/2., y0-l/2.), (x0+t/2., y0+l/2.), layer))  
  
def allignment_mark_offset(x0, y0, size, cw, cross_layer, dark_layer, light_layer, gap, offset):
        cross(x0-1.*size, y0, cw+1., size+cw, cross_layer)
        cross(x0+1.*size, y0, cw+1., size+cw, cross_layer)
        box_grid(x0-1.*size, y0+offset, (size-2.*gap+cw)/2., 2*gap, 2, light_layer)
        box_grid(x0+1.*size, y0+offset, (size-2.*gap+cw)/2., 2.*gap, 2, dark_layer)
        poly_cell.add(gdspy.Text(str(offset), 100, (x0, y0+size/3.+100), layer = cross_layer))
        allignment_mark_box(x0, y0, cross_layer, w = 5.*size, t = size/2.)
        allignment_mark_box(x0, y0, dark_layer, w = 5.*size, t = size/2.)
        allignment_mark_box(x0, y0, light_layer, w = 5.*size, t = size/2.)
        
def fine_alignment_mark(x0, y0, gap, layer1, layer2, n2 = 6):
        nodd = 2*n2 + 1
        t = gap
        l = 5.*gap
        s_offset = l
        lag = (t+gap)/2.
        #vertical marks 
        for i in range(nodd):
                print i
                if  i != n2 - 1 and i != n2 + 1 and i != n2:
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 - l, y0 + s_offset + i*(t + gap) + t ), layer1))
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 - l, y0 + s_offset + i*(t + gap) + t ), layer1).rotate(2.*numpy.pi-numpy.pi/2., center = [x0, y0]))
                        
                        
                        if i>n2 :
                                poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap) - lag -t/4.), (x0 - l, y0 + s_offset + i*(t + gap) + t - lag + t/4.), layer2))
                                poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap) - lag-t/4.), (x0 - l, y0 + s_offset + i*(t + gap) + t - lag +t/4.), layer2).rotate(2.*numpy.pi-numpy.pi/2., center = [x0, y0]))
                        if i<n2:
                                poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap) + lag - t/4.), (x0 - l, y0 + s_offset + i*(t + gap) + lag +t + t/4.), layer2))
                                poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap) + lag - t/4.), (x0 - l, y0 + s_offset + i*(t + gap) + lag +t + t/4.), layer2).rotate(2.*numpy.pi-numpy.pi/2., center = [x0, y0]))
                        
                if i == n2:
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 - 2.*l, y0 + s_offset + i*(t + gap) + t), layer1))
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 - 2.*l, y0 + s_offset + i*(t + gap) + t), layer2))
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 + 2.*l, y0 + s_offset + i*(t + gap) + t), layer1).rotate(2.*numpy.pi-numpy.pi/2., center = [x0, y0]))
                        poly_cell.add(gdspy.Rectangle((x0 , y0 + s_offset + i*(t + gap)), (x0 + 2.*l, y0 + s_offset + i*(t + gap) + t), layer2).rotate(2.*numpy.pi-numpy.pi/2., center = [x0, y0]))
                

fine_alignment_mark(-36500, -700, 40, 3, 1)
fine_alignment_mark(35500, -700, 40, 3, 1)                        


centy = -(6.*7500.+5000.)/2.
centx = -(10.*3750.+1850.)/2.
cantilever_grid(centx, centy, 11, 7, 3750, 7500)
allignment_mark_box(-40000, 0, 1, w = 2000)
allignment_mark_box(40000, 0, 1, w = 2000)

allignment_mark_box(-40000, 0, 2, w = 2000)
allignment_mark_box(40000, 0, 2, w = 2000)

allignment_mark_box(-40000, 0, 3, w = 2000)
allignment_mark_box(40000, 0, 3, w = 2000)


#text = gdspy.Text('Sample text', 2000, (-10, -100))
#poly_cell.add(text)

allignment_mark_offset(42500, 0, 100, 10, 3, 2, 1, 4, 2)
allignment_mark_offset(-42500, 0, 100, 10, 3, 2, 1, 4, 2)
allignment_mark_offset(37500, 0, 100, 10, 3, 2, 1, 4, -2)
allignment_mark_offset(-37500, 0, 100, 10, 3, 2, 1, 4, -2)
#connection
#trench_pattern(0, 0, 25, 100, 100, 5000, 500, 10, 1)
## ------------------------------------------------------------------ ##
##	OUTPUT															  ##
## ------------------------------------------------------------------ ##

## Output the layout to a GDSII file (default to all created cells).
## Set the units we used to micrometers and the precision to nanometers.



## ------------------------------------------------------------------ ##
##	IMPORT															  ##
## ------------------------------------------------------------------ ##


gdsii1 = gdspy.GdsImport('Cross_marks_27.gds', unit = 1e-6, layers = {1:3, 2:2, 50:50}, rename = {'base_for_clear_mark':'base_for_clear_mark1', 'clear_layer_mark':'clear_layer_mark1'})

gdsii1.extract('4in_karlsuss_backside')
gdsii1.extract('base_for_clear_mark1')
gdsii1.extract('clear_layer_mark1')
gdsii1.extract('4in_wafer')

poly_cell.add(gdspy.CellReference('4in_karlsuss_backside', (0, 0)))
poly_cell.add(gdspy.CellReference('base_for_clear_mark1', (40000-245+30+500, 750)))
poly_cell.add(gdspy.CellReference('clear_layer_mark1', (40000-245+30+500, 750)))

poly_cell.add(gdspy.CellReference('base_for_clear_mark1', (-40000-245+30-500, 750)))
poly_cell.add(gdspy.CellReference('clear_layer_mark1', (-40000-245+30-500, 750)))

gdsii3 = gdspy.GdsImport('Cross_marks_27.gds', unit = 1e-6, layers = {1:3, 2:1, 50:50}, rename = {'base_for_clear_mark':'base_for_clear_mark3', 'clear_layer_mark':'clear_layer_mark3', '4in_wafer':'pos', '4in_karlsuss_backside':'shitty'})


gdsii3.extract('base_for_clear_mark3')
gdsii3.extract('clear_layer_mark3')



poly_cell.add(gdspy.CellReference('base_for_clear_mark3', (40000-245+30-500, 750)))
poly_cell.add(gdspy.CellReference('clear_layer_mark3', (40000-245+30-500, 750)))

poly_cell.add(gdspy.CellReference('base_for_clear_mark3', (-40000-245+30+500, 750)))
poly_cell.add(gdspy.CellReference('clear_layer_mark3', (-40000-245+30+500, 750)))


gdsii2 = gdspy.GdsImport('ks_MA6_ff_4in.gds', layers = {1:3, 2:50})
gdsii2.extract('Cell0')
poly_cell.add(gdspy.CellReference('Cell0', (0, 0)))


#ref_cell.add(gdspy.CellReference('ff', (0, 0)))

gdspy.gds_print('attractorv3.gds', unit=1.0e-6, precision=1.0e-9)
## ------------------------------------------------------------------ ##
##	VIEWER															  ##
## ------------------------------------------------------------------ ##


## View the layout using a GUI.  Full description of the controls can
## be found in the online help at http://gdspy.sourceforge.net/
gdspy.LayoutViewer()
