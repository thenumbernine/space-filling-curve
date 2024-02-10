#!/usr/bin/env luajit
require 'ext'
local bit = require 'bit'

-- TODO I competley forgot what I wanted to do with this project
-- something about average distances between points of space-filling curves
-- maybe to see which ones were better for coherency


--[==[
-- n = hilbert iteration
-- i = integer [0,2^(2n))
-- returns x,y in [0,2^n)^2
local function hilberti(i,n)
	assert(n >= 0)
	assert(i >= 0 and i < bit.lshift(1,2*n))
	if n == 0 then return 0, 0 end
	if n == 1 then
		local bx = bit.rshift(bit.band(i,2),1)
		local by = bit.bxor(bit.band(i,1),bx)
		return by, bx
	end
	local sx,sy = hilberti(bit.band(i, bit.lshift(1,2*(n-1))-1), n-1)
end

local t = table()
local n = 2
for i=0,bit.lshift(1,2*n)-1 do
	local x,y = hilberti(i,n)
	x=x+1 y=y+1
	x,y = y,x
	y = bit.lshift(1,n)-y+1
	t[y] = t[y] or table()
	t[y][x] = i
end
print(t:map(function(l) return l:concat'\t' end):concat'\n')


-- t is in 0-1
-- n is the hilbert curve iteration
-- returns (x,y) in [0,1]^2
local function hilbert(t,n)
	local nSq = n * n
	t = t * nSq
	local i = math.floor(t)
	local frac = t - i
	if i == nSq then
		i = nSq - 1
		frac = 1
	end

	local x1,y1 = hilberti(i,n)
	local x2,y2 = hilberti(i+1,n)
	return x1 + (x2 - x1) * frac, y1 + (y2 - y1) * frac
end
--]==]
-- want signedness
--local vec2l = require 'vec-ffi.create_vec2'{ctype='int64_t', suffix='l'}
local vec2d = require 'vec-ffi.vec2d'

local Curve = class()


-- https://en.wikipedia.org/wiki/Z-order_curve
local ZCurve = Curve:subclass()

ZCurve.name = 'Z curve' 

function ZCurve:build(iter)
	local pts = table()
	local n = bit.lshift(1, bit.lshift(iter, 1))
	for i=0,n-1 do
		local x = bit.bor(
			bit.band(2^0, i),
			bit.band(2^1, bit.rshift(i, 1)),
			bit.band(2^2, bit.rshift(i, 2)),
			bit.band(2^3, bit.rshift(i, 3)),
			bit.band(2^4, bit.rshift(i, 4)),
			bit.band(2^5, bit.rshift(i, 5)),
			bit.band(2^6, bit.rshift(i, 6)),
			bit.band(2^7, bit.rshift(i, 7))
		)
		local y = bit.bor(
			bit.band(2^0, bit.rshift(i, 1)),
			bit.band(2^1, bit.rshift(i, 2)),
			bit.band(2^2, bit.rshift(i, 3)),
			bit.band(2^3, bit.rshift(i, 4)),
			bit.band(2^4, bit.rshift(i, 5)),
			bit.band(2^5, bit.rshift(i, 6)),
			bit.band(2^6, bit.rshift(i, 7)),
			bit.band(2^7, bit.rshift(i, 8))
		)
		pts:insert(vec2d(x, y))
	end
	return pts
end


-- https://en.wikipedia.org/wiki/L-system
local RewriteCurve = Curve:subclass()
RewriteCurve.turnAngle = 90
RewriteCurve.lengths = {f = 1, g = 1}
function RewriteCurve:build(iter)
	local s = self.axiom
	for i=0,iter-1 do
		s = s:gsub('.', self.rules)
	end
	local pts = table()
	local th = math.rad(self.turnAngle)
	local r = vec2d(math.cos(th), math.sin(th))
	local n = vec2d(1,0)
	local p = vec2d(0,0)
	pts:insert(p:clone())
	for i=1,#s do
		local c = s:sub(i,i)
		local l = self.lengths[c]
		if l then
			p = p + l * n
		end
		if c == '+' then
			--n.x, n.y = -n.y, n.x	-- 90' opt
			n.x, n.y =
				r.x * n.x - r.y * n.y,
				r.x * n.y + r.y * n.x
		elseif c == '-' then
			--n.x, n.y = n.y, -n.x
			n.x, n.y =
				r.x * n.x + r.y * n.y,
				r.x * n.y - r.y * n.x
		end
		pts:insert(p:clone())
	end
	return pts
end


-- https://en.wikipedia.org/wiki/Hilbert_curve
local HilbertCurve = RewriteCurve:subclass()
HilbertCurve.name = 'Hilbert curve' 
HilbertCurve.axiom = 'a'
HilbertCurve.rules = {
	a = '+bf-afa-fb+',
	b = '-af+bfb+fa-',
}


-- https://en.wikipedia.org/wiki/Moore_curve
local MooreCurve = RewriteCurve:subclass()
MooreCurve.name = 'Moore curve' 
MooreCurve.axiom = 'lfl+f+lfl'
MooreCurve.rules = {
	l = '-rf+lfl+fr-',
	r = '+lf-rfr-fl+',
}
-- one iteration off
function MooreCurve:build(iter)
	return MooreCurve.super.build(self, iter-1)
end


-- https://en.wikipedia.org/wiki/L-system
local KochCurve = RewriteCurve:subclass()
KochCurve.name = 'Koch curve'
KochCurve.axiom = 'f'
KochCurve.rules = {
	f = 'f+f-f-f+f',
}


-- https://en.wikipedia.org/wiki/L-system
local SierpinskiTriangle = RewriteCurve:subclass()
SierpinskiTriangle.name = 'Sierpinski triangle'
SierpinskiTriangle.turnAngle = 120
SierpinskiTriangle.axiom = 'f-g-g'
SierpinskiTriangle.rules = {
	f = 'f-g+f+g-f',
	g = 'gg',
}

-- https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve
local SierpinskiCurve = RewriteCurve:subclass()
SierpinskiCurve.name = 'Sierpinski curve'
SierpinskiCurve.turnAngle = 45
SierpinskiCurve.axiom = 'f--xf--f--xf'
SierpinskiCurve.rules = {
	x = 'xf+g+xf--f--xf+g+x',
}

-- https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve
local SierpinskiSquareCurve = RewriteCurve:subclass()
SierpinskiSquareCurve.name = 'Sierpinski square curve'
SierpinskiSquareCurve.axiom = 'f+xf+f+xf'
SierpinskiSquareCurve.rules = {
	x = 'xf-f+f-xf+f+xf-f+f-x',
}

-- https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve
local SierpinskiArrowHeadCurve = RewriteCurve:subclass()
SierpinskiArrowHeadCurve.name = 'Sierpinski arrowhead curve'
SierpinskiArrowHeadCurve.turnAngle = 60
SierpinskiArrowHeadCurve.axiom = 'xf'
SierpinskiArrowHeadCurve.rules = {
	x = 'yf+xf+y',
	y = 'xf-yf-x',
}

local curves = table{
	ZCurve(),
	HilbertCurve(),
	MooreCurve(),
	KochCurve(),
	SierpinskiTriangle(),
	SierpinskiCurve(),
	SierpinskiSquareCurve(),
	SierpinskiArrowHeadCurve(),
	-- PeanoCurve() -- https://en.wikipedia.org/wiki/Peano_curve

}
local curveNames = curves:mapi(function(c) return c.name end)

local gl = require 'gl'
local ig = require 'imgui'
local App = require 'imguiapp.withorbit'()

App.title = 'space filling curves'

function App:initGL(...)
	App.super.initGL(self, ...)
	self.view.ortho = true
	self.view.orthoSize = 1
	self.curve = 1
	self.iter = 1
	self.tex = require 'gl.hsvtex2d'(256, nil, true)
		:unbind()
	self:rebuild()
end

function App:rebuild()
	self.iter = math.clamp(self.iter, 1, 8)
	self.path = curves[self.curve]:build(self.iter)
end

function App:update(...)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	local s = 1 / tonumber(bit.lshift(1, self.iter))
	gl.glTranslatef(-.5, -.5, 0)
	gl.glScalef(s, s, 1)
	self.tex
		:enable()
		:bind()
	gl.glBegin(gl.GL_LINE_STRIP)
	for i,p in ipairs(self.path) do
		gl.glTexCoord1f((i-.5)/#self.path)
		gl.glVertex2f(p:unpack())
	end
	gl.glEnd()
	self.tex
		:unbind()
		:disable()
	App.super.update(self, ...)
end

function App:updateGUI()
	if ig.luatableCombo('curve', self, 'curve', curveNames)
	or ig.luatableInputInt('iter', self, 'iter')
	then
		self:rebuild()
	end
end

return App():run()
