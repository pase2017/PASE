# 所有 .c 文件
file(GLOB_RECURSE lib_src *.c)
# 所有 .h 文件
file(GLOB_RECURSE lib_header *.h)
# 拷贝 .h 文件至指定路径
file(COPY ${lib_header} DESTINATION ${pase_include_output_path})
# 添加库文件
add_library(${pase_library_hypre} STATIC ${lib_src})