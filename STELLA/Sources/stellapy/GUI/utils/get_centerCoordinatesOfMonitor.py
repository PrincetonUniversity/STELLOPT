
import screeninfo
def get_centerCoordinatesOfMonitor(root, window, width=None, heigth=None):
    ''' When using dual monitors, this will give the center of the monitor where the root is located. '''

    # Get the monitors
    monitors = screeninfo.get_monitors()
    
    # Get the location of the root window
    x = root.winfo_x()
    y = root.winfo_y()

    # Get the coordiantes of the center of the first screen in case our loop does not find any matches
    center_x = monitors[0].x + (monitors[0].width/2)
    center_y = monitors[0].y + (monitors[0].height/2)    
        
    # Decide which monitor corresponds to (x,y)
    for m in reversed(monitors):
        if m.x <= x <= m.width + m.x and m.y <= y <= m.height + m.y:
            # Get the coordinates of the center of this monitor
            center_x = m.x + (m.width/2)
            center_y = m.y + (m.height/2)
    
    # Substract the dimensions of our new toplevel window
    # To ensure the new toplevel is centered on the screen
    x = center_x - window.winfo_reqwidth()/2    if width==None else  center_x - width/2
    y = center_y - window.winfo_reqheight()/2   if width==None else  center_y - heigth/2
    return x, y

