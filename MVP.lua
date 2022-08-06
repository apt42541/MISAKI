function WorldToScreen(s, w, mat)
    local vec = Vec3TransformCoord(w, mat)
    s[1] = (1.0 + vec[1]) * 1280.0 / 2.0
    s[2] = (1.0 - vec[2]) * 720.0 / 2.0
    s[3] = vec[3]

    if vec[3] < 1.0 then
       return true
    end

    return false
 end

local mvpmat = MatrixRead(GetAddress("ac_client.exe+17DFD0"))
local pos = Vec3()
local screenPos = Vec2()
pos[1] = 113.1239471
pos[2] = 84.41746521
pos[3] = 0.0

if WorldToScreen(screenPos, pos, mvpmat) then
print(string.format("true: %f, x: %f, y: %f", screenPos[3], screenPos[1], screenPos[2]))
else
print(string.format("false: %f", screenPos[3]))
end
