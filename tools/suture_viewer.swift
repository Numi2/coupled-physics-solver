import AppKit
import SceneKit

struct Frame: Decodable {let tissue:[[Float]];let lips:[[[Float]]]?;let thread:[[[Float]]];let gap:Double;let pull:Double;let tension:Double;let time:Double}
struct Replay: Decodable {let triangles:[Int32];let frames:[Frame]}

final class Viewer: NSObject, NSApplicationDelegate, NSWindowDelegate {
    let replay:Replay; let output:String?
    let window=NSWindow(contentRect:NSRect(x:0,y:0,width:1000,height:760),styleMask:[.titled,.closable,.resizable,.miniaturizable],backing:.buffered,defer:false)
    let view=SCNView();let scene=SCNScene();let tissue=SCNNode();let strand=SCNNode()
    let label=NSTextField(labelWithString:"");let slider=NSSlider();let button=NSButton(title:"Pause",target:nil,action:nil)
    var timer:Timer?;var frame=0;var playing=true;var xray=true
    init(replay:Replay,output:String?) {self.replay=replay;self.output=output;super.init()}
    func applicationDidFinishLaunching(_ notification:Notification) {
        window.title="Skin Suture — Metal replay";window.delegate=self
        let root=NSView();window.contentView=root
        view.scene=scene;view.backgroundColor=NSColor(calibratedRed:0.94,green:0.95,blue:0.97,alpha:1);view.allowsCameraControl=true;view.autoenablesDefaultLighting=true;view.antialiasingMode = .multisampling4X
        view.translatesAutoresizingMaskIntoConstraints=false;root.addSubview(view)
        let title=NSTextField(labelWithString:"Synthetic skin · five prethreaded stitches · prescribed tail pull")
        let subtitle=NSTextField(labelWithString:"Elastic thread slides through ideal eyelets. Skin positions come from the native FEM solve.")
        title.font = .systemFont(ofSize:18,weight:.semibold);subtitle.textColor = .secondaryLabelColor
        button.target=self;button.action=#selector(toggle)
        slider.minValue=0;slider.maxValue=Double(replay.frames.count-1);slider.target=self;slider.action=#selector(scrub)
        let opacity=NSButton(checkboxWithTitle:"See thread inside",target:self,action:#selector(toggleOpacity));opacity.state = .on
        let controls=NSStackView(views:[button,slider,opacity]);controls.orientation = .horizontal;slider.widthAnchor.constraint(greaterThanOrEqualToConstant:400).isActive=true
        let footer=NSStackView(views:[title,subtitle,label,controls]);footer.orientation = .vertical;footer.alignment = .leading;footer.spacing=10;footer.translatesAutoresizingMaskIntoConstraints=false;root.addSubview(footer)
        NSLayoutConstraint.activate([view.topAnchor.constraint(equalTo:root.topAnchor),view.leadingAnchor.constraint(equalTo:root.leadingAnchor),view.trailingAnchor.constraint(equalTo:root.trailingAnchor),view.bottomAnchor.constraint(equalTo:footer.topAnchor,constant:-12),footer.leadingAnchor.constraint(equalTo:root.leadingAnchor,constant:20),footer.trailingAnchor.constraint(equalTo:root.trailingAnchor,constant:-20),footer.bottomAnchor.constraint(equalTo:root.bottomAnchor,constant:-20)])
        scene.rootNode.addChildNode(tissue);scene.rootNode.addChildNode(strand)
        let camera=SCNNode();camera.camera=SCNCamera();camera.camera!.usesOrthographicProjection=true;camera.camera!.orthographicScale=20;camera.camera!.zNear=0.1;camera.camera!.zFar=200
        camera.position=SCNVector3(12,-28,30);camera.look(at:SCNVector3(0,0,0));scene.rootNode.addChildNode(camera);view.pointOfView=camera
        show(0);window.center();window.makeKeyAndOrderFront(nil);NSApp.activate(ignoringOtherApps:true)
        if let output=output {
            for i in replay.frames.indices {show(i)}
            DispatchQueue.main.asyncAfter(deadline:.now()+1) {
                let image=self.view.snapshot();let data=image.tiffRepresentation!;let rep=NSBitmapImageRep(data:data)!
                try! rep.representation(using:.png,properties:[:])!.write(to:URL(fileURLWithPath:output))
                print("native_suture_viewer=ok frames=\(self.replay.frames.count) final_frame=\(self.frame) screenshot=\(output)")
                NSApp.terminate(nil)
            }
        } else {timer=Timer.scheduledTimer(withTimeInterval:1.0/30,repeats:true) {_ in if self.playing {self.show((self.frame+1)%self.replay.frames.count)}}}
    }
    // Display tubes follow exported coordinates; radius is a visualization
    // choice, not a calibrated thread diameter or an added physics constraint.
    func tube(_ a:SCNVector3,_ b:SCNVector3,radius:CGFloat,color:NSColor)->SCNNode? {
        let dx=b.x-a.x,dy=b.y-a.y,dz=b.z-a.z;let length=sqrt(dx*dx+dy*dy+dz*dz)
        if length<1e-6 {return nil}
        let geometry=SCNCylinder(radius:radius,height:CGFloat(length));geometry.radialSegmentCount=8
        geometry.firstMaterial?.diffuse.contents=color;geometry.firstMaterial?.lightingModel = .blinn
        let node=SCNNode(geometry:geometry);node.position=SCNVector3((a.x+b.x)/2,(a.y+b.y)/2,(a.z+b.z)/2)
        let axis=sqrt(dx*dx+dz*dz)
        node.rotation=axis>1e-6 ? SCNVector4(dz/axis,0,-dx/axis,acos(max(-1,min(1,dy/length)))):SCNVector4(1,0,0,dy<0 ? Float.pi:0)
        return node
    }
    func show(_ i:Int) {
        frame=i;let f=replay.frames[i]
        let vertices=f.tissue.map{SCNVector3($0[0],$0[1],$0[2])}
        let indices=replay.triangles;let element=SCNGeometryElement(indices:indices,primitiveType:.triangles)
        let geometry=SCNGeometry(sources:[SCNGeometrySource(vertices:vertices)],elements:[element])
        let skin=SCNMaterial();skin.diffuse.contents=NSColor(calibratedRed:0.79,green:0.48,blue:0.47,alpha:1);skin.isDoubleSided=true;skin.lightingModel = .blinn;skin.roughness.contents=0.85;skin.transparency=xray ? 0.40:1;skin.writesToDepthBuffer = !xray
        geometry.materials=[skin];tissue.geometry=geometry
        strand.childNodes.forEach{$0.removeFromParentNode()}
        for points in f.thread {
            let vv=points.map{SCNVector3($0[0],$0[1],$0[2])}
            for j in 0..<(vv.count-1) {
                if let node=tube(vv[j],vv[j+1],radius:0.035,color:.systemBlue) {strand.addChildNode(node)}
            }
            for j in stride(from:0,to:vv.count,by:4) {let dot=SCNSphere(radius:0.055);dot.firstMaterial?.diffuse.contents=NSColor.systemBlue;let node=SCNNode(geometry:dot);node.position=vv[j];strand.addChildNode(node)}
        }
        for points in f.lips ?? [] {
            let vv=points.map{SCNVector3($0[0],$0[1],$0[2])}
            if vv.count>1 {for j in 0..<(vv.count-1) {
                if let node=tube(vv[j],vv[j+1],radius:0.012,color:NSColor(calibratedRed:0.35,green:0.10,blue:0.12,alpha:1)) {strand.addChildNode(node)}
            }}
        }
        slider.doubleValue=Double(i)
        label.stringValue=String(format:"Frame %d/%d · %.1f ms · mean gap %.3f mm · pull %.3f mm · peak stitch tension %.2f mN",i+1,replay.frames.count,f.time*1000,f.gap,f.pull,f.tension)
    }
    @objc func toggleOpacity(){xray.toggle();show(frame)}
    @objc func toggle(){playing.toggle();button.title=playing ? "Pause":"Play"}
    @objc func scrub(){playing=false;button.title="Play";show(Int(slider.doubleValue.rounded()))}
    func windowWillClose(_ notification:Notification){timer?.invalidate();NSApp.terminate(nil)}
}
let args=CommandLine.arguments
if args.count<2 {fputs("usage: Suture replay.json [--snapshot output.png]\n",stderr);exit(2)}
do {
    let replay=try JSONDecoder().decode(Replay.self,from:Data(contentsOf:URL(fileURLWithPath:args[1])))
    guard !replay.frames.isEmpty else {throw NSError(domain:"empty replay",code:1)}
    let delegate=Viewer(replay:replay,output:args.count==4 && args[2]=="--snapshot" ? args[3]:nil)
    let app=NSApplication.shared;app.setActivationPolicy(.regular);app.delegate=delegate;app.run()
} catch {fputs("suture_viewer=failed \(error)\n",stderr);exit(1)}
