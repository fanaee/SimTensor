About SimTensor
-----------------

SimTensor software has been developed at the Laboratory of Artificial Intelligence in INESC TEC research institute in the context of European FP7 Project "MAESTRA". The goal of SimTensor project is to provide a multi-platform, open-source software for generating artificial tensor-structured data (CP/PARAFAC and Tucker) for reproducible research on tensor factorization algorithms. SimTensor is a stand-alone application based on MATALB (You do not need a MATLAB license to run it). It provides a wide range of facilities for generating tensor data with various configurations. SimTensor comes with a user-friendly graphical user interface, which enables the user to generate tensors with complicated settings in an extremely easy and fast way. It also has this facility to export generated data to universal formats such as CSV and HDF5, which can be imported via a wide range of programming languages (C, C++, Java, R, Fortran, MATLAB, Perl, Python, and many more). The most innovative part of SimTensor is that you can generate temporal tensors with periodic waves, seasonal effects and streaming structure. You also can apply constraints such as non-negativity and different kinds of sparsity to your data. SimTensor also gives you the option to simulate different kinds of change-points and inject various types of anomalies in your data.


How to run
----------
1. Download Tensor toolbox: http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
2. Extract tensor_toolbox_2.6.zip into root directory
3. Run UI.m to start SimTensor

How to cite
-----------

If you used the SimTensor in your work in any way, please cite the following paper:

Fanaee-T, H., & Gama, J. (2016). SimTensor: A synthetic tensor data generator. arXiv preprint arXiv:1406.3496.

```
@article{fanaee2016simtensor,
	title={ SimTensor: A synthetic tensor data generator},
	author={Fanaee-T, Hadi and Gama, Joao},
	journal={arXiv preprint arXiv:1612.03772},
	year={2016}
}
```

Compile note
------------
If you want to compile to binary do not forget to add `*.txt` files inside of the `GUIt/` directory


Team
----
Hadi Fanaee-T, University of Porto, Portugal
Joao Gama, University of Porto, Portugal

Credit: Andrzej Cichocki, Anh-Huy Phan, Arun Tejasvi Chaganty, Brett W. Bader, Changwei Hu, Changyou Chen, Christos Faloutsos, Claus A. Andersson, Daniel M. Dunlavy, Dimitri Nion, Eric C. Chi, Evrim Acar, Giorgio Tomasi, Guoxu Zhou, Haesun Park, Ian Davidson, Jackson Mayo, James Bailey, Jimeng Sun, Jingu Kim, Laurent Sorber, Lawrence Carin, Lieven De Lathauwer, Liqing Zhang, Marc Van Barel, Masatoshi Yoshikawa, Matthew Harding, Nguyen Xuan Vinh, Nicholas D. Sidiropoulos, Nico Vervliet, Otto Debals, Papalexakis, Evangelos E., Percy Liang, Petr Tichavsky, Piyush Rai, Prateek Jain, Qibin Zhao, Rafal Zdunekb, Rasmus Bro, Rong Pan, Sammy Hansen, Sewoong Oh, Shun-ichi Amari, Shuo Zhou, Tamara G. Kolda, Tatsuya Yokotaa, Tomoharu Iwata, Volodymyr Kuleshov, Wotao Yin, Xiaomin Fang, Xingyu Wang , Yangyang Xu, Yasuko Matsubara, Yasushi Sakurai, Yu Zhang , Yukihiko Yamashitac, Yunlong He, Yunzhe Jia


Feedback and bug report
-----------------------
Feedbacks are very welcome. Thank you!
email: [support@simtensor.org](mailto: support@simtensor.org)
