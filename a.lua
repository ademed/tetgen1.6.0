-- Lua script.
p=tetview:new()
<<<<<<< HEAD
p:load_mesh("barout")
=======
p:load_mesh("c.1")
>>>>>>> 464193039a1791db0d7e708b0d14f27a6226f5bd
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
