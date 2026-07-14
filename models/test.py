import json

candi_pts = []

with open(r"../data/test_candidate_points005.json", "r") as f:
    candi_pts = json.load(f)

pts = []
for i in candi_pts:
    for j in i:
        if j != []:
            pts.append(j)

def get_x_boundary(points, index):
    '''
    Takes a set of points and returns the domain of the set
    :param points: list of (x,y) pairs
    :return: the maximum and minimum x values of the set
    '''
    xmax = points[0][index]
    for i in range(len(points)):
        if points[i][index] > xmax:
            xmax = points[i][index]
    xmin = points[0][index]
    for i in range(len(points)):
        if points[i][index] < xmin:
            xmin = points[i][index]
    return xmin, xmax

print("x bounds")
print(get_x_boundary(pts,0))
print("y bounds")
print(get_x_boundary(pts,1))
print("z bounds")
print(get_x_boundary(pts,2))


if [-150,50,22.5] in pts:
    print("num found")
