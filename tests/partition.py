
import math

def partition2x2(ulx, uly, lrx, lry):
    n = 1
    partitions = []
    remaining = [(ulx, uly, lrx, lry)]

    while remaining:
        rect = remaining.pop()
        width = rect[2] - rect[0] + 1 
        height = rect[3] - rect[1] + 1

        # Use whichever dimension is smallest as base size of square
        if width < height:
            n = math.floor(math.pow(2, (math.floor(math.log(width, 2)))))
            n = int(n)
        else:
            n = math.floor((math.pow(2, (math.floor(math.log(height, 2))))))
            n = int(n)
        # n is now equal to the width of the square and a power of two
        # Append it to the paritions list
        partitions.append((rect[0], rect[1], rect[0]+n-1, rect[1]+n-1))

        if (rect[0]+n-1) != rect[2]:
            # handle leftover right side
            remaining.append((rect[0]+n, rect[1], rect[2], rect[3]))

        if (rect[1]+n-1) != rect[3]:
            # handle leftover bottom
            remaining.append((rect[0], rect[1]+n, rect[0]+n-1, rect[3]))

    return partitions

def get_partition_area(partition):
    width = partition[2] - partition[0] + 1
    height = partition[3] - partition[1] + 1

    return width * height
    
def survey_areas(partitions):
    areas = {}

    for p in partitions:
        a = get_partition_area(p)
        if a in areas:
            areas[a] += 1
        else:
            areas[a] = 1

    return areas

def merge_rectangle_list(rectangles):
    while len(rectangles) != 1:
        r1 = rectangles.pop(0)
        for i,v in enumerate(rectangles):
            new = merge_rectangles(r1, v)
            if new:
                rectangles[i] = new
                continue 

def merge_rectangles(rect1, rect2):
    """
    If the two arguments share a side, merge them and return a tuple
    with the upper-left and lower-right corners of the new rectangle.
    If the two arguments do not share a side, return None.
    """
    if rect2 == None:
        return None
    
    if rectangle_above(rect1, rect2):
        return (rect1[0], rect1[1], rect2[2], rect2[3])
    elif rectangle_left(rect1, rect2):
        return (rect1[0], rect1[1], rect2[2], rect2[3])
    elif rectangle_above(rect2, rect1):
        return (rect2[0], rect2[1], rect1[2], rect1[3])
    elif rectangle_left(rect2, rect1):
        return (rect2[0], rect2[1], rect1[2], rect1[3])
    else:
        return None

def rectangle_above(rect1, rect2):
    """ Check whether the two arguments share a side and the shared
    side the is first's bottom and the second's top, i.e. rect1 is
    above rect2.
    """
    if rect1[3] == rect2[1]-1:
        # rect1 is exactly 1 above rect2
        if rect1[0] == rect2[0]:
                # rect1 ulx == rect2 ulx
                if rect1[2] == rect2[2]:
                    # rect1 lrx == rect2 lrx
                    return True


    return False

def rectangle_left(rect1, rect2):
    """ Check whether the two arguments share a side and the shared
    side the is first's right and the second's left, i.e. rect1 is
    left of rect2.
    """
    if rect1[2] == rect2[0] - 1:
        if rect1[1] == rect2[1]:
                if rect1[3] == rect2[3]:
                    return True


    return False

