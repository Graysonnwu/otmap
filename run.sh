# clear


./build/transport_points -in_src data/source.png -in_trg data/yyqx.png -res 800 -point_light -light_pos 0 0 0.6 -focal 0 -thick 0.05 -o pixlens0322-brick40-s1-f06.obj
# ./build/transport_points -in_src data/source.png -in_trg data/zkg.png -res 300 -point_light -light_pos 0 0 1 -trg_scale 4 4 1 -focal 3 -thick 0.1 -o test_thickness.obj
# ./build/transport_points -in_src data/source.png -in_trg data/paojie2.png -res 1200 -focal 3 -o paojie2.obj
# ./build/transport_points -in_src data/source-head.png -in_trg data/head2.png -res 500 -focal 3 -o head.obj
# ./build/transport_points -in_src data/source.png -in_trg data/zkg.png -res 500 -light_pos _1 0 1.854777 -trg_pos 1 0 _1 -trg_rot 0 _90 0 -o zkg-2.obj
# ./build/transport_points -in_src data/source.png -in_trg data/zkg.png -res 500 -light_pos 1 0 1.854777 -trg_pos _1 0 _1 -trg_rot 0 90 0 -o zkg.obj
# ./build/transport_points -in_src data/source.png -in_trg data/zkg.png -res 400 -light_pos 0 0 1 -trg_pos 0 0 _2 -trg_rot 0 0 0 -o zkg.obj

# ./build/transport_points -in_src data/source.png -in_trg data/leaf.png -res 500 -light_pos _1 0 1.854777 -trg_pos 2 0 _2 -trg_rot 0 _90 0 -o leaf4.obj # 45°折射
# ./build/transport_points -in_src data/source.png -in_trg data/leaf.png -res 500 -trg_pos 4 0 _4 -trg_scale 1.5 1.5 1 -o leaf2.obj -reflect # 45°反射
# ./build/transport_points -in_src data/source.png -in_trg data/leaf.png -res 500 -light_pos 0 0 _1 -trg_pos 0 0 _3 -trg_rot 0 0 0 -o leaf3.obj -reflect # 垂直反射