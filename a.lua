-- Lua script.
p=tetview:new()
p:load_mesh("C:/Users/ADEYEMI AHMED/Downloads/tetgen1.6.0/tetgen1.6.0/pmdc.1.ele")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
