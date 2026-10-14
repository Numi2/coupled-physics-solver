#include "../matter/src/matter_kernels.metal"
using namespace numi_matter_metal;
kernel void suture_ccd_contract(device float4* result [[buffer(0)]], uint i [[thread_position_in_grid]]) {
    float toi=0,sep=INFINITY;float3 normal,point,bary;bool certified=false,hit=false;
    float3 a(-.01,-.01,0),b(.01,-.01,0),c(0,.01,0);
    if(i==0) { // Common translation must not manufacture relative speed.
        float3 move(.1,.02,.03),v(0,0,.000020);
        hit=vertexTriangleCCD(v,v+move,a,a+move,b,b+move,c,c+move,.000010,toi,normal,bary,point,sep,certified);
    } else if(i==1) { // A true vertex/triangle crossing.
        hit=vertexTriangleCCD(float3(0,0,.001),float3(0,0,-.001),a,a,b,b,c,c,.000010,toi,normal,bary,point,sep,certified);
    } else if(i==2) { // Common translation of parallel edges.
        float3 move(.1,.02,.03),p(-.001,0,0),q(.001,0,0),r(-.001,.000020,0),s(.001,.000020,0);float2 param;
        hit=edgeEdgeCCD(p,p+move,q,q+move,r,r+move,s,s+move,.000010,toi,normal,param,point,sep,certified);
    } else if(i==3) { // A true crossing of perpendicular segments.
        float2 param;
        hit=edgeEdgeCCD(float3(-.001,0,0),float3(-.001,0,0),float3(.001,0,0),float3(.001,0,0),float3(0,-.001,.001),float3(0,-.001,-.001),float3(0,.001,.001),float3(0,.001,-.001),.000010,toi,normal,param,point,sep,certified);
    }
    if(i>=4 && i<7) {
        float scale=i==4?1.0f:(i==5?.001f:.00001f);float3 g0,g1;
        float value=ipcFeatureMollifierDifferential(float3(scale,0,0),float3(0,scale,0),g0,g1);
        result[i]=float4(value,length(g0),length(g1),1);return;
    }
    if(i==7) {
        float3 bary;
        float3 closest=closestPointOnTriangle(float3(.00025,.00025,.00001),float3(0,0,0),float3(.001,0,0),float3(0,.00075,0),bary);
        result[i]=float4(closest,bary.y);return;
    }
    if(i==8) {
        float s,t;float3 p,q;
        closestPointsOnSegments(float3(-.00025,0,0),float3(.00025,0,0),float3(0,-.00025,.00001),float3(0,.00025,.00001),s,t,p,q);
        result[i]=float4(s,t,length(p-q),1);return;
    }
    result[i]=float4(hit?1:0,certified?1:0,toi,sep);
}
