function Matrix()
   local m = {}

   for i = 0, 15 do
      local v = 0.0

      if (i % 5) == 0 then
         v = 1.0
      end

      table.insert(m, v)
   end

   return m
end

function MatrixMultiply(m1, m2)
   local m = Matrix()

   for i = 0, 3 do
      for j = 0, 3 do
         m[i*4+j+1] = m1[i*4+1] * m2[j+1] + m1[i*4+2] * m2[4+j+1] + m1[i*4+3] * m2[8+j+1] + m1[i*4+4] * m2[12+j+1]
      end
   end

   return m
end

function MatrixRead(p)
   local m = Matrix()

   for i = 0, 15 do
      m[i+1] = ReadFloat(p + i * 0x4)
   end

   return m
end

function MatrixInverse(m)
    local det = 0.0
    local t = Vec3()
    local v = Matrix()
    local out = Matrix()

    t[1] = m[11] * m[16] - m[12] * m[15]
    t[2] = m[7] * m[16] - m[8] * m[15]
    t[3] = m[7] * m[12] - m[8] * m[11]
    v[1] = m[6] * t[1] - m[10] * t[2] + m[14] * t[3]
    v[5] = -m[5] * t[1] + m[9] * t[2] - m[13] * t[3]

    t[1] = m[5] * m[10] - m[9] * m[6]
    t[2] = m[5] * m[14] - m[13] * m[6]
    t[3] = m[9] * m[14] - m[13] * m[10]
    v[9] = m[16] * t[1] - m[12] * t[2] + m[8] * t[3]
    v[13] = -m[15] * t[1] + m[11] * t[2] - m[7] * t[3]

    det = m[1] * v[1] + m[2] * v[5] + m[3] * v[9] + m[4] * v[13]

    if det == 0.0 then
        return Matrix()
    end

    t[1] = m[11] * m[16] - m[12] * m[15]
    t[2] = m[3] * m[16] - m[4] * m[15]
    t[3] = m[3] * m[12] - m[4] * m[11]
    v[2] = -m[2] * t[1] + m[10] * t[2] - m[14] * t[3]
    v[6] = m[1] * t[1] - m[9] * t[2] + m[13] * t[3]

    t[1] = m[1] * m[10] - m[9] * m[2]
    t[2] = m[13] * m[2] - m[1] * m[14]
    t[3] = m[9] * m[14] - m[13] * m[10]
    v[10] = -m[16] * t[1] - m[12] * t[2] - m[4] * t[3]
    v[14] = m[15] * t[1] + m[11] * t[2] + m[3] * t[3]

    t[1] = m[7] * m[16] - m[8] * m[15]
    t[2] = m[3] * m[16] - m[4] * m[15]
    t[3] = m[3] * m[8] - m[4] * m[7]
    v[3] = m[2] * t[1] - m[6] * t[2] + m[14] * t[3]
    v[7] = -m[1] * t[1] + m[5] * t[2] - m[13] * t[3]

    t[1] = m[1] * m[6] - m[5] * m[2]
    t[2] = m[13] * m[2] - m[1] * m[14]
    t[3] = m[5] * m[14] - m[13] * m[6]
    v[11] = m[16] * t[1] + m[8] * t[2] + m[4] * t[3]
    v[15] = -m[15] * t[1] - m[7] * t[2] - m[3] * t[3]

    t[1] = m[7] * m[12] - m[8] * m[11]
    t[2] = m[3] * m[12] - m[4] * m[11]
    t[3] = m[3] * m[8] - m[4] * m[7]
    v[4] = -m[2] * t[1] + m[6] * t[2] - m[10] * t[3]
    v[8] = m[1] * t[1] - m[5] * t[2] + m[9] * t[3]

    v[12] = -m[1] * (m[6] * m[12] - m[8] * m[10]) + m[5] * (m[2] * m[12] - m[4] * m[10]) - m[9] * (m[2] * m[8] - m[4] * m[6])
    v[16] = m[1] * (m[6] * m[11] - m[7] * m[10]) - m[5] * (m[2] * m[11] - m[3] * m[10]) + m[9] * (m[2] * m[7] - m[3] * m[6])

    det = 1.0 / det

    for i = 0, 4 - 1 do
         for j = 0, 4 - 1 do
            out[i*4+j+1] = v[4 * i + j + 1] * det
	      end
    end

    return out
end

function MatrixInvert(M)
   local s = {}
	local c = {}
   local T = Matrix()

	s[1] = M[1]*M[6] - M[5]*M[2]
	s[2] = M[1]*M[7] - M[5]*M[3]
	s[3] = M[1]*M[8] - M[5]*M[4]
	s[4] = M[2]*M[7] - M[6]*M[3]
	s[5] = M[2]*M[8] - M[6]*M[4]
	s[6] = M[3]*M[8] - M[7]*M[4]

	c[1] = M[9]*M[14] - M[13]*M[10]
	c[2] = M[9]*M[15] - M[13]*M[11]
	c[3] = M[9]*M[16] - M[13]*M[12]
	c[4] = M[10]*M[15] - M[14]*M[11]
	c[5] = M[10]*M[16] - M[14]*M[12]
	c[6] = M[11]*M[16] - M[15]*M[12]

	-- Assumes it is invertible
	local idet = 1.0 / ( s[1]*c[6]-s[2]*c[5]+s[3]*c[4]+s[4]*c[3]-s[5]*c[2]+s[6]*c[1] )

	T[1] = ( M[6] * c[6] - M[7] * c[5] + M[8] * c[4]) * idet
	T[2] = (-M[2] * c[6] + M[3] * c[5] - M[4] * c[4]) * idet
	T[3] = ( M[14] * s[6] - M[15] * s[5] + M[16] * s[4]) * idet
	T[4] = (-M[10] * s[6] + M[11] * s[5] - M[12] * s[4]) * idet

	T[5] = (-M[5] * c[6] + M[7] * c[3] - M[8] * c[2]) * idet
	T[6] = ( M[1] * c[6] - M[3] * c[3] + M[4] * c[2]) * idet
	T[7] = (-M[13] * s[6] + M[15] * s[3] - M[16] * s[2]) * idet
	T[8] = ( M[9] * s[6] - M[11] * s[3] + M[12] * s[2]) * idet

	T[9] = ( M[5] * c[5] - M[6] * c[3] + M[8] * c[1]) * idet
	T[10] = (-M[1] * c[5] + M[2] * c[3] - M[4] * c[1]) * idet
	T[11] = ( M[13] * s[5] - M[14] * s[3] + M[16] * s[1]) * idet
	T[12] = (-M[9] * s[5] + M[10] * s[3] - M[12] * s[1]) * idet

	T[13] = (-M[5] * c[4] + M[6] * c[2] - M[7] * c[1]) * idet
	T[14] = ( M[1] * c[4] - M[2] * c[2] + M[3] * c[1]) * idet
	T[15] = (-M[13] * s[4] + M[14] * s[2] - M[15] * s[1]) * idet
	T[16] = ( M[9] * s[4] - M[10] * s[2] + M[11] * s[1]) * idet

   return T
end

function MatrixTranspose(m)
   local out = Matrix()

   for i = 0, 3 do
       for j = 0, 3 do
          out[i*4+j+1] = m[j*4+i+1]
       end
   end

   return out
end

function MatrixRotationYawPitchRoll(yaw, pitch, roll)
   local sroll, croll, spitch, cpitch, syaw, cyaw
   local out = Matrix()

   sroll = math.sin(roll)
   croll = math.cos(roll)
   spitch = math.sin(pitch)
   cpitch = math.cos(pitch)
   syaw = math.sin(yaw)
   cyaw = math.cos(yaw)

   out[1] = sroll * spitch * syaw + croll * cyaw
   out[2] = sroll * cpitch
   out[3] = sroll * spitch * cyaw - croll * syaw
   out[4] = 0.0
   out[5] = croll * spitch * syaw - sroll * cyaw
   out[6] = croll * cpitch
   out[7] = croll * spitch * cyaw + sroll * syaw
   out[8] = 0.0
   out[9] = cpitch * syaw
   out[10] = -spitch
   out[11] = cpitch * cyaw
   out[12] = 0.0
   out[13] = 0.0
   out[14] = 0.0
   out[15] = 0.0
   out[16] = 1.0

   return out
end

function MatrixLookAtLH(eye, at, up)
    local right   = Vec3()
    local upn     = Vec3()
    local forward = Vec3()
    local out     = Matrix()

    forward       = Vec3Subtract(at, eye)
    forward       = Vec3Normalize(forward)
    right         = Vec3Cross(up, forward)
    upn           = Vec3Cross(forward, right)
    right         = Vec3Normalize(right)
    upn           = Vec3Normalize(upn)

    out[1]        = right[1]
    out[5]        = right[2]
    out[9]        = right[3]
    out[13]       = -Vec3Dot(right, eye)
    out[2]        = upn[1]
    out[6]        = upn[2]
    out[10]       = upn[3]
    out[14]       = -Vec3Dot(upn, eye)
    out[3]        = forward[1]
    out[7]        = forward[2]
    out[11]       = forward[3]
    out[15]       = -Vec3Dot(forward, eye)
    out[4]        = 0.0
    out[8]        = 0.0
    out[12]       = 0.0
    out[16]       = 1.0

    return out
end

function MatrixLookAtRH(eye, at, up)
    local right   = Vec3()
    local upn     = Vec3()
    local forward = Vec3()
    local out     = Matrix()

    vec           = Vec3Subtract(at, eye)
    vec           = Vec3Normalize(vec)
    right         = Vec3Cross(up, vec)
    upn           = Vec3Cross(vec, right)
    right         = Vec3Normalize(right)
    upn           = Vec3Normalize(upn)

    out[1]        = -right[1]
    out[5]        = -right[2]
    out[9]        = -right[3]
    out[13]       = Vec3Dot(right, eye)
    out[2]        = upn[1]
    out[6]        = upn[2]
    out[10]       = upn[3]
    out[14]       = -Vec3Dot(upn, eye)
    out[3]        = -vec[1]
    out[7]        = -vec[2]
    out[11]       = -vec[3]
    out[15]       = Vec3Dot(vec, eye)
    out[4]        = 0.0
    out[8]        = 0.0
    out[12]       = 0.0
    out[16]       = 1.0

    return out
end

function MatrixTranslation(m, v)
    m[13] = v[1]
    m[14] = v[2]
    m[15] = v[3]

    return m
end

function MatrixGetTranslation(m)
    local out = Vec3()

    out[1] = m[13]
    out[2] = m[14]
    out[3] = m[15]

    return out
end

function MatrixGetRotation(m)
   local out = Matrix()

   for i = 0, 2 do
      for l = 0, 2 do
          out[i*4+l+1] = m[i*4+l+1]
      end
   end

   return out
end

function MatrixPrint(m)
   local s = ''

   for i = 0, 15 do
      if i > 0 then
         if (i % 4) == 0 then
            s = s .. '\n'
         else
            s = s .. ' '
         end
      end

      s = s.. string.format('%f', m[i + 1])
   end

   print(s)
end









function Vec2()
   local v = {}

   for i = 0, 1 do
       table.insert(v, 0.0)
   end

   return v
end

function Vec2Read(p)
    local v = Vec3()

    for i = 0, 1 do
       v[i + 1] = ReadFloat(p + i * 0x4)
    end

    return v
end

function Vec2Length(v)
    return math.sqrt(Vec2LengthSq(v))
end

function Vec2LengthSq(v)
    return v[1] * v[1] + v[2] * v[2]
end











function Vec3()
   local v = {}

   for i = 0, 2 do
      table.insert(v, 0.0)
   end

   return v
end

function Vec3Read(p)
   local v = Vec3()

   for i = 0, 2 do
      v[i + 1] = ReadFloat(p + i * 0x4)
   end

   return v
end

function Vec3Negative(v)
   v[1] = -v[1]
   v[2] = -v[2]
   v[3] = -v[3]

   return v
end

function Vec3Cross(v1, v2)
   local out = Vec3()

   out[1] = v1[2] * v2[3] - v1[3] * v2[2]
   out[2] = v1[3] * v2[1] - v1[1] * v2[3]
   out[3] = v1[1] * v2[2] - v1[2] * v2[1]

   return out
end

function Vec3Dot(v1, v2)
   return v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
end

function Vec3TransformCoord(v, m)
   local out = Vec3()
   local norm = m[4] * v[1] + m[8] * v[2] + m[12] * v[3] + m[16]
   out[1] = (m[1] * v[1] + m[5] * v[2] + m[9] * v[3] + m[13]) / norm
   out[2] = (m[2] * v[1] + m[6] * v[2] + m[10] * v[3] + m[14]) / norm
   out[3] = (m[3] * v[1] + m[7] * v[2] + m[11] * v[3] + m[15]) / norm
   return out
end

function Vec3Length(v)
   return math.sqrt(Vec3LengthSq(v))
end

function Vec3LengthSq(v)
   return v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
end

function Vec3Normalize(v)
   local out = Vec3()
   local norm = Vec3Length(v)

   if not norm then
      out[1] = 0.0
      out[2] = 0.0
      out[3] = 0.0
   else
      out[1] = v[1] / norm
      out[2] = v[2] / norm
      out[3] = v[3] / norm
   end
end

function Vec3Subtract(v1, v2)
    local out = Vec3()

    out[1] = v1[1] - v2[1]
    out[2] = v1[2] - v2[2]
    out[3] = v1[3] - v2[3]

    return out
end

function Vec3Transform(v, m)
    local out = Vec4()

    out[1] = m[1] * v[1] + m[5] * v[2] + m[9] * v[3] + m[13]
    out[2] = m[2] * v[1] + m[6] * v[2] + m[10] * v[3] + m[14]
    out[3] = m[3] * v[1] + m[7] * v[2] + m[11] * v[3] + m[15]
    out[4] = m[4] * v[1] + m[8] * v[2] + m[12] * v[3] + m[16]

    return out
end

function Vec3TransformNormal(v, m)
    local out = Vec3()

    out[1] = m[1] * v[1] + m[5] * v[2] + m[9] * v[3]
    out[2] = m[2] * v[1] + m[6] * v[2] + m[10] * v[3]
    out[3] = m[3] * v[1] + m[7] * v[2] + m[11] * v[3]

    return out
end

function Vec3Print(v)
   print(string.format("%f %f %f", v[1], v[2], v[3]))
end





function Vec4()
   local v = {}

   for i = 0, 4 - 1 do
      table.insert(v, 0.0)
   end

   return v
end

function Vec4Transform(v, m)
    local out = Vec4()

    out[1] = m[1] * v[1] + m[5] * v[2] + m[9] * v[3] + m[13] * v[4]
    out[2] = m[2] * v[1] + m[6] * v[2] + m[10] * v[3] + m[14] * v[4]
    out[3] = m[3] * v[1] + m[7] * v[2] + m[11] * v[3] + m[15] * v[4]
    out[4] = m[4] * v[1] + m[8] * v[2] + m[12] * v[3] + m[16] * v[4]

    return out
end






function WorldToScreen(s, w, mat)
    local vec = Vec3TransformCoord(w, mat)
    s[2] = (1.0 + vec[2]) * 1920.0 / 2.0
    s[3] = (1.0 - vec[3]) * 1080.0 / 2.0

    if vec[4] < 1.0 then
       return true
    end

    return false
 end

 function MakeProjection(fov, aspect)
    local out = Vec2()
    local tanFov = math.tan(math.rad(fov * .5))

    out[1] = 1.0 / (tanFov * aspect)
    out[2] = 1.0 / tanFov

    return out
 end

 function InverseProjection(v)
    local out = Vec2()

    out[2] = v[2] / v[1]
    out[1] = math.deg(2.0 * math.atan(1.0 / v[1] / out[2]))

    return out
 end
