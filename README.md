# NS-2D Solver
此repository包括偏微分离散，超大线性方程组迭代等核心模块的代码和测试，不是完整代码，完整代码未开源。
只是给行业从业者和爱好者提供参考和交流。

### 效果展示
1. `3D`
   - [The Curl Distribution in Rotor Wake](http://v.youku.com/v_show/id_XMTY0NzM1MDQyMA==.html)
2. `2D`
   - [Pressure Distribution of Dynamic Stall in Reverse Flow](http://v.youku.com/v_show/id_XMTYxOTU0MzQ5Mg==.html)</br>
   - [Kármán Vortex Street after Aerodynamic Trailing Edge](http://v.youku.com/v_show/id_XMTYxOTU0MzgyNA==.html)</br>
   - [Structured Grid Generation around NACA0012](http://v.youku.com/v_show/id_XMTYxOTU0NDUzMg==.html)</br>
   - [Unstructured Grid Generation outside Structured Grid](http://v.youku.com/v_show/id_XMTYxOTU0NDc2OA==.html)

### keys:
- steady: 只计算静态
- unsteady: 可以计算静态和动态
- Convective flux: JST ; Roe
- Temporal Discretization: Runge Kutta; Dual-Time+LU-SGS
- Turbulent Model: SA; k- $\\omega$ SST
- Moving-grid
