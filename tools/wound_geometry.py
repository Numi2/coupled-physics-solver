"""Independent double-precision wound-face intersection observer."""
import math

def sub(a,b):return [x-y for x,y in zip(a,b)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def segment_triangle(p,q,a,b,c):
    edge1=sub(b,a);edge2=sub(c,a);d=sub(q,p);h=cross(d,edge2);det=dot(edge1,h)
    if abs(det)<1e-24:return False
    inv=1/det;s=sub(p,a);u=inv*dot(s,h)
    if u<0 or u>1:return False
    z=cross(s,edge1);v=inv*dot(d,z)
    if v<0 or u+v>1:return False
    t=inv*dot(edge2,z);return 0<t<1

def winding(a,b,c,p):return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
def triangle_overlap(a,b):
    # Non-coplanar intersection by six segment/triangle tests.
    if any(segment_triangle(a[i],a[(i+1)%3],*b) or segment_triangle(b[i],b[(i+1)%3],*a) for i in range(3)):return True
    n=cross(sub(a[1],a[0]),sub(a[2],a[0]));length=math.sqrt(dot(n,n))
    if length==0:return False
    # A narrow coplanar path closes the usual segment/plane blind spot.
    if any(abs(dot(n,sub(v,a[0])))/length>1e-12 for v in b):return False
    drop=max(range(3),key=lambda d:abs(n[d]));project=lambda v:[v[d] for d in range(3) if d!=drop]
    aa=[project(v) for v in a];bb=[project(v) for v in b]
    def inside(p,t):
        s=[winding(t[i],t[(i+1)%3],p) for i in range(3)]
        return all(v>1e-20 for v in s) or all(v<-1e-20 for v in s)
    if any(inside(p,bb) for p in aa) or any(inside(p,aa) for p in bb):return True
    for i in range(3):
        for j in range(3):
            p,q=aa[i],aa[(i+1)%3];r,s=bb[j],bb[(j+1)%3]
            if winding(p,q,r)*winding(p,q,s)<0 and winding(r,s,p)*winding(r,s,q)<0:return True
    return False

def check_wound(nodes,faces):
    rest=nodes[min(nodes)];lo=[t for t in faces if all(abs(rest[n][1]+.0003)<1e-8 for n in t)];hi=[t for t in faces if all(abs(rest[n][1]-.0003)<1e-8 for n in t)]
    assert lo and hi,'missing actual incision faces'
    tests=0
    for step,positions in nodes.items():
        for aa in lo:
            a=[positions[n] for n in aa];alo=[min(v[d] for v in a) for d in range(3)];ahi=[max(v[d] for v in a) for d in range(3)]
            for bb in hi:
                if set(aa)&set(bb):continue
                b=[positions[n] for n in bb]
                if any(max(v[d] for v in b)<alo[d]-1e-12 or min(v[d] for v in b)>ahi[d]+1e-12 for d in range(3)):continue
                tests+=1
                assert not triangle_overlap(a,b),f'wound faces intersect at frame {step}: {aa}, {bb}'
    return {'wound_surface_intersections':0,'wound_surface_candidate_tests':tests,'wound_faces_per_side':[len(lo),len(hi)]}
