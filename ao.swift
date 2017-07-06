#if os(Linux)
import Glibc
#else
import Darwin
#endif

let NAOSamples = 8
let NSubSamples = 2
let width = 256
let height = 256

struct Plane {
    let p : Vector
    let n : Vector
    
    init(p: Vector, n: Vector) {
        self.p = p
        self.n = n
    }
    
    func rayIntersect(isect : Isect, ray : Ray) -> Isect
    {
        let d = -self.p.dot(v : self.n);
        let v = ray.dir.dot(v : self.n);
        
        if fabs(v) < 1.0e-17 {
            return isect
        }
        
        let t = -(ray.org.dot(v : self.n) + d) / v;
        
        if (t > 0.0) && (t < isect.t) {
            let p = Vector(x: ray.org.x + ray.dir.x * t,
                           y: ray.org.y + ray.dir.y * t,
                           z: ray.org.z + ray.dir.z * t)
            
            return Isect(t: t, p: p, n: self.n, hit: true)
        }
        
        return isect
    }
}

struct Sphere {
    let center : Vector
    let radius : Double
    
    init(center: Vector, radius: Double) {
        self.center = center
        self.radius = radius
    }
    
    func rayIntersect(isect : Isect, ray : Ray) -> Isect {
        let rs = Vector(x: ray.org.x - self.center.x,
                        y: ray.org.y - self.center.y,
                        z: ray.org.z - self.center.z)
        
        
        let B = rs.dot(v : ray.dir)
        let C = rs.dot(v : rs) - self.radius * self.radius;
        let D = B * B - C
        
        if D > 0.0 {
            let t = -B - sqrt(D)
            
            if (t > 0.0) && (t < isect.t) {
                let p = Vector(x: ray.org.x + ray.dir.x * t,
                               y: ray.org.y + ray.dir.y * t,
                               z: ray.org.z + ray.dir.z * t)
                let n = Vector(x: p.x - self.center.x,
                               y: p.y - self.center.y,
                               z: p.z - self.center.z)
                return Isect(t: t, p: p, n: n.normalize(), hit: true)
            }
        }
        
        return isect
    }
}

struct Isect {
    let t : Double
    let p : Vector
    let n : Vector
    let hit : Bool
    
    init(t: Double, p: Vector, n: Vector, hit: Bool) {
        self.t = t
        self.p = p
        self.n = n
        self.hit = hit
    }
}

struct Ray {
    let org : Vector
    let dir : Vector
    
    init(org: Vector, dir: Vector) {
        self.org = org
        self.dir = dir
    }
}

struct Vector {
    let x : Double
    let y : Double
    let z : Double
    
    init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
    
    static let zero = Vector(x : 0.0, y: 0.0, z: 0.0)
    
    func dot(v : Vector) -> Double {
        return self.x * v.x + self.y * v.y + self.z * v.z
    }
    
    func cross(v : Vector) -> Vector {
        return Vector(x: self.y * v.z - self.z * v.y,
                      y: self.z * v.x - self.x * v.z,
                      z: self.x * v.y - self.y * v.x)
    }
    
    func length() -> Double {
        return sqrt(self.dot(v: self))
    }
    
    func normalize() -> Vector {
        let length = self.length()
        
        if fabs(length) > 1.0e-17 {
            return Vector(x: self.x / length,
                          y: self.y / length,
                          z: self.z / length)
        }
        
        return self
    }
}

extension Vector {
    typealias Basis = [Vector]
    
    func orthoBasis() -> Basis {
        var basis = Basis(repeating: Vector.zero, count : 3)
        basis[2] = self;
        
        switch (self.x, self.y, self.z) {
        case (let x, _, _) where x < 0.6 && x > -0.6:
            basis[1] = Vector(x : 1.0, y: 0.0, z: 0.0)
        case (_, let y, _) where y < 0.6 && y > -0.6:
            basis[1] = Vector(x : 0.0, y: 1.0, z: 0.0)
        case (_, _, let z) where z < 0.6 && z > -0.6:
            basis[1] = Vector(x : 0.0, y: 0.0, z: 1.0)
        default:
            basis[1] = Vector(x : 1.0, y: 0.0, z: 0.0)
        }
        
        basis[0] = basis[1].cross(v: basis[2]).normalize()
        basis[1] = basis[2].cross(v: basis[0]).normalize()
        
        return basis
    }
}

extension Double {
    func clamp() -> UInt8 {
        switch Int32(self * 255.5) {
        case let i where i < 0:
            return 0
        case let i where i > 255:
            return 255
        default:
            return UInt8(self * 255.5)
        }
    }   
}

let spheres = [
    Sphere(center: Vector(x: -2.0, y: 0.0, z: -3.5),
           radius: 0.5),
    Sphere(center: Vector(x: -0.5, y: 0.0, z: -3.0),
           radius: 0.5),
    Sphere(center: Vector(x: 1.0, y: 0.0, z: -2.2),
           radius: 0.5),
    Sphere(center: Vector(x: -2.0, y: 0.0, z: 0.0),
           radius: 0.5)
]
let plane = Plane(p: Vector(x: 0.0, y: -0.5, z: 0.0),
                  n: Vector(x: 0.0, y: 1.0 , z: 0.0))

func ambientOcclusion(isect : Isect) -> Vector
{
    let ntheta : Int = NAOSamples
    let nphi : Int = NAOSamples
    let eps : Double = 0.0001
    
    let p = Vector(x: isect.p.x + eps * isect.n.x,
                   y: isect.p.y + eps * isect.n.y,
                   z: isect.p.z + eps * isect.n.z)
    
    let basis = isect.n.orthoBasis()
    var occlusion : Double = 0.0;
    
    for _ in 0..<ntheta {
        for _ in 0..<nphi {
            let theta = sqrt(drand48());
            let phi = 2.0 * Double.pi * drand48();
            
            let x = cos(phi) * theta;
            let y = sin(phi) * theta;
            let z = sqrt(1.0 - theta * theta);
            
            // local . global
            let rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
            let ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
            let rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;
            
            let ray = Ray(org: p, dir: Vector(x: rx, y: ry, z: rz))
            var occIsect = Isect(t: 1.0e+17, p: Vector.zero, n: Vector.zero, hit: false)
            
            occIsect = spheres[0].rayIntersect(isect: occIsect, ray: ray)
            occIsect = spheres[1].rayIntersect(isect: occIsect, ray: ray)
            occIsect = spheres[2].rayIntersect(isect: occIsect, ray: ray)
            occIsect = spheres[3].rayIntersect(isect: occIsect, ray: ray)
            occIsect = plane.rayIntersect(isect: occIsect, ray: ray)
            
            if occIsect.hit {
                occlusion += 1.0
            }
        }
    }
    occlusion = (Double(ntheta * nphi) - occlusion) / Double(ntheta * nphi)

    return Vector(x: occlusion, y: occlusion, z: occlusion)
}

func render(w : Int, h : Int, nsubsamples : Int) -> [UInt8]
{
    var img = [UInt8](repeating: 0, count: w * h * 3)
    var fimg = [Double](repeating: 0, count: w * h * 3)
    
    for y in 0..<h {
        for x in 0..<w {
            for v in 0..<nsubsamples {
                for u in 0..<nsubsamples     {
                    let w2 = Double(w) / 2.0
                    let h2 = Double(h) / 2.0
                    let px = (Double(x) + (Double(u) / Double(nsubsamples)) - w2) / w2
                    let py = -(Double(y) + (Double(v) / Double(nsubsamples)) - h2) / h2
                    
                    let ray = Ray(org: Vector(x: 0.0, y: 0.0, z: 0.0),
                                  dir: Vector(x: px, y: py, z: -1.0).normalize())
                    var isect = Isect(t: 1.0e+17, p: Vector.zero, n: Vector.zero, hit: false)
                    isect = spheres[0].rayIntersect(isect: isect, ray: ray)
                    isect = spheres[1].rayIntersect(isect: isect, ray: ray)
                    isect = spheres[2].rayIntersect(isect: isect, ray: ray)
                    isect = spheres[3].rayIntersect(isect: isect, ray: ray)
                    isect = plane.rayIntersect(isect: isect, ray: ray)
                    
                    if isect.hit {
                        let col = ambientOcclusion(isect: isect);
                        
                        fimg[3 * (y * w + x) + 0] += col.x
                        fimg[3 * (y * w + x) + 1] += col.y
                        fimg[3 * (y * w + x) + 2] += col.z
                    }
                    
                }
            }
            
            fimg[3 * (y * w + x) + 0] /= Double(nsubsamples * nsubsamples)
            fimg[3 * (y * w + x) + 1] /= Double(nsubsamples * nsubsamples)
            fimg[3 * (y * w + x) + 2] /= Double(nsubsamples * nsubsamples)
            
            img[3 * (y * w + x) + 0] = (fimg[3 * (y * w + x) + 0]).clamp()
            img[3 * (y * w + x) + 1] = (fimg[3 * (y * w + x) + 1]).clamp()
            img[3 * (y * w + x) + 2] = (fimg[3 * (y * w + x) + 2]).clamp()
        }
    }
    
    return img
}

func saveppm(fname : String, w : Int, h : Int, img : [UInt8])
{
    guard let fp = fopen(fname, "wb") else {
        print("Open Error \(strerror(errno))")
        return
    }
    
    let _ = withVaList([]) {
        vfprintf(fp, "P6\n", $0);
    }
    let _ = withVaList([w, h]) {
        vfprintf(fp, "%d %d\n", $0);        
    }
    let _ = withVaList([]) {
        vfprintf(fp, "255\n", $0);        
    }
    fwrite(img, w * h * 3, 1, fp);
    fclose(fp);
}

let img = render(w: width, h: height, nsubsamples: NSubSamples);
saveppm(fname: "ao.ppm", w: width, h: height, img: img);
