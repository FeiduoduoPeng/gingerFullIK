#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Geometry>

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#define mySign(x) ((x>0)-(x<0))

typedef Eigen::Matrix<double, 10, 1> Vector10d;
typedef Eigen::Matrix<double, 10, 10> Matrix10d;

class gingerFullIK {

public:
    gingerFullIK(){
        M23_ << cos(TILT), -sin(TILT), 0, 0.25743421,
                sin(TILT), cos(TILT), 0, -0.167643,
                0, 0, 1, -0.01532,
                0, 0, 0, 1 ;
        M45_ << cos(-TILT), -sin(-TILT), 0, -0.22791,
                sin(-TILT), cos(-TILT), 0, 0.03067462,
                0, 0, 1, 0,
                0, 0, 0, 1 ;
        M67_ << 1, 0, 0, -0.209557780,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;
        cstr << -1.54, 1.54, //waist
                -0.50, 0.20,
                -0.38, 0.38,
                -0.60, 3.14, //shoudler
                -0.38, 1.54,
                -1.54, 1.54,
                -0.20, 2.10, //elbow
                -1.54, 1.54, //wrist
                -0.38, 0.38,
                -0.60, 0.38;
    }

private:
    const double PI = 3.14159265359;
    const double TILT = 20.0/180.0*PI;
    const int MAX_ITER = 2000;
    double dt = 0.1;
    double condThresh = 200.0;
    bool optFlag = false;
    Eigen::Matrix4d M01, M12, M23, M34, M45, M56, M67, M78, M89, M910;
    Eigen::Matrix4d M010;
    Eigen::Matrix4d M10_9, M10_8, M10_7, M10_6, M10_5, M10_4, M10_3, M10_2, M10_1, M10_0;
    Eigen::Matrix4d M23_, M45_, M67_;
public:
    Eigen::Matrix<double, 10, 2> cstr; //constraint

private:
    void setM(const Vector10d &angs){
        double th0 = angs(0), th1 = angs(1), th2 = angs(2), th3 = angs(3), th4 = angs(4),
               th5 = angs(5), th6 = angs(6), th7 = angs(7), th8 = angs(8), th9 = angs(9);

        M01 << 1, 0, 0, 0,
               0, cos(th0), -sin(th0), 0,
               0, sin(th0), cos(th0), 0,
               0, 0, 0, 1;
        M12 << cos(th1), 0, sin(th1), 0,
               0, 1, 0, 0,
               -sin(th1), 0, cos(th1), 0,
               0, 0, 0, 1;
        M23 << cos(th2), -sin(th2), 0, 0.105,
               sin(th2), cos(th2), 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;
        M34 << cos(th3), 0,  sin(th3), 0,
               0, 1, 0, 0,
               -sin(th3), 0, cos(th3), 0,
               0, 0, 0, 1 ;
        M45 << cos(th4), -sin(th4), 0, 0,
               sin(th4), cos(th4), 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;
        M56 << 1, 0, 0, 0,
               0, cos(th5), -sin(th5), 0,
               0, sin(th5), cos(th5), 0,
               0, 0, 0, 1;
        M67 << cos(th6), 0, sin(th6),  0,
               0, 1, 0, 0,
               -sin(th6), 0, cos(th6),  0,
               0, 0, 0, 1;
        M78 << 1, 0, 0, 0,
               0, cos(th7), -sin(th7), 0,
               0, sin(th7), cos(th7), 0,
               0, 0, 0, 1;
        M89 << cos(th8), 0,  sin(th8), 0,
               0, 1, 0, 0,
               -sin(th8), 0, cos(th8), 0,
               0, 0, 0, 1;
        M910 << cos(th9), -sin(th9), 0, 0,
                sin(th9), cos(th9), 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

        M23 = M23*M23_;
        M45 = M45*M45_;
        M67 = M67*M67_;

        M010 = M01*M12*M23*M34*M45*M56*M67*M78*M89*M910;

        M10_9 = M910.inverse();
        M10_8 = M10_9 * (M89.inverse());
        M10_7 = M10_8 * (M78.inverse());
        M10_6 = M10_7 * (M67.inverse());
        M10_5 = M10_6 * (M56.inverse());
        M10_4 = M10_5 * (M45.inverse());
        M10_3 = M10_4 * (M34.inverse());
        M10_2 = M10_3 * (M23.inverse());
        M10_1 = M10_2 * (M12.inverse());
        M10_0 = M10_1 * (M01.inverse());
    }

    Vector10d utheta2Theta(const Vector10d &utheta){
        return Vector10d(
            (cstr(0,1)-cstr(0,0)) / PI * atan(utheta(0)) + (cstr(0,1)+cstr(0,0)) / 2.0,
            (cstr(1,1)-cstr(1,0)) / PI * atan(utheta(1)) + (cstr(1,1)+cstr(1,0)) / 2.0,
            (cstr(2,1)-cstr(2,0)) / PI * atan(utheta(2)) + (cstr(2,1)+cstr(2,0)) / 2.0,
            (cstr(3,1)-cstr(3,0)) / PI * atan(utheta(3)) + (cstr(3,1)+cstr(3,0)) / 2.0,
            (cstr(4,1)-cstr(4,0)) / PI * atan(utheta(4)) + (cstr(4,1)+cstr(4,0)) / 2.0,
            (cstr(5,1)-cstr(5,0)) / PI * atan(utheta(5)) + (cstr(5,1)+cstr(5,0)) / 2.0,
            (cstr(6,1)-cstr(6,0)) / PI * atan(utheta(6)) + (cstr(6,1)+cstr(6,0)) / 2.0,
            (cstr(7,1)-cstr(7,0)) / PI * atan(utheta(7)) + (cstr(7,1)+cstr(7,0)) / 2.0,
            (cstr(8,1)-cstr(8,0)) / PI * atan(utheta(8)) + (cstr(8,1)+cstr(8,0)) / 2.0,
            (cstr(9,1)-cstr(9,0)) / PI * atan(utheta(9)) + (cstr(9,1)+cstr(9,0)) / 2.0
        );
    }

    Vector10d theta2UTheta(const Vector10d &theta){
        return Vector10d (
            tan( PI/2 * (2*theta(0)-cstr(0,1)-cstr(0,0)) / (cstr(0,1)-cstr(0,0)) ),
            tan( PI/2 * (2*theta(1)-cstr(1,1)-cstr(1,0)) / (cstr(1,1)-cstr(1,0)) ),
            tan( PI/2 * (2*theta(2)-cstr(2,1)-cstr(2,0)) / (cstr(2,1)-cstr(2,0)) ),
            tan( PI/2 * (2*theta(3)-cstr(3,1)-cstr(3,0)) / (cstr(3,1)-cstr(3,0)) ),
            tan( PI/2 * (2*theta(4)-cstr(4,1)-cstr(4,0)) / (cstr(4,1)-cstr(4,0)) ),
            tan( PI/2 * (2*theta(5)-cstr(5,1)-cstr(5,0)) / (cstr(5,1)-cstr(5,0)) ),
            tan( PI/2 * (2*theta(6)-cstr(6,1)-cstr(6,0)) / (cstr(6,1)-cstr(6,0)) ),
            tan( PI/2 * (2*theta(7)-cstr(7,1)-cstr(7,0)) / (cstr(7,1)-cstr(7,0)) ),
            tan( PI/2 * (2*theta(8)-cstr(8,1)-cstr(8,0)) / (cstr(8,1)-cstr(8,0)) ),
            tan( PI/2 * (2*theta(9)-cstr(9,1)-cstr(9,0)) / (cstr(9,1)-cstr(9,0)) )
        );
    }

    Matrix10d uthetaJ(Vector10d utheta){
        // theta = (U-L)/pi/*atan(utheta) + (U+L)/2. So derive is: (U-L)/ pi / (utheta^2+1)
        Vector10d diag(
            (cstr(0,1)-cstr(0,0)) / (utheta(0)*utheta(0)+1.0) / PI,
            (cstr(1,1)-cstr(1,0)) / (utheta(1)*utheta(1)+1.0) / PI,
            (cstr(2,1)-cstr(2,0)) / (utheta(2)*utheta(2)+1.0) / PI,
            (cstr(3,1)-cstr(3,0)) / (utheta(3)*utheta(3)+1.0) / PI,
            (cstr(4,1)-cstr(4,0)) / (utheta(4)*utheta(4)+1.0) / PI,
            (cstr(5,1)-cstr(5,0)) / (utheta(5)*utheta(5)+1.0) / PI,
            (cstr(6,1)-cstr(6,0)) / (utheta(6)*utheta(6)+1.0) / PI,
            (cstr(7,1)-cstr(7,0)) / (utheta(7)*utheta(7)+1.0) / PI,
            (cstr(8,1)-cstr(8,0)) / (utheta(8)*utheta(8)+1.0) / PI,
            (cstr(9,1)-cstr(9,0)) / (utheta(9)*utheta(9)+1.0) / PI
        );
        Matrix10d uthetaJ = diag.asDiagonal(); //uthetaJ is positive-define

        return uthetaJ;
    }

    Eigen::Matrix<double, 6, 10> gingerJ(const Vector10d &angs){
        setM(angs);
        Eigen::Matrix<double, 6, 10> Jfull_10 = Eigen::Matrix<double, 6, 10>::Zero();

        Jfull_10.block<3,1>(0,0) = M10_0.block<3,1>(0,3).cross( M10_0.block<3,1>(0,0) );
        Jfull_10.block<3,1>(0,1) = M10_1.block<3,1>(0,3).cross( M10_1.block<3,1>(0,1) );
        Jfull_10.block<3,1>(0,2) = M10_2.block<3,1>(0,3).cross( M10_2.block<3,1>(0,2) );
        Jfull_10.block<3,1>(0,3) = M10_3.block<3,1>(0,3).cross( M10_3.block<3,1>(0,1) );
        Jfull_10.block<3,1>(0,4) = M10_4.block<3,1>(0,3).cross( M10_4.block<3,1>(0,2) );
        Jfull_10.block<3,1>(0,5) = M10_5.block<3,1>(0,3).cross( M10_5.block<3,1>(0,0) );
        Jfull_10.block<3,1>(0,6) = M10_6.block<3,1>(0,3).cross( M10_6.block<3,1>(0,1) );
        Jfull_10.block<3,1>(0,7) = M10_7.block<3,1>(0,3).cross( M10_7.block<3,1>(0,0) );
        Jfull_10.block<3,1>(0,8) = M10_8.block<3,1>(0,3).cross( M10_8.block<3,1>(0,1) );
        Jfull_10.block<3,1>(0,9) = M10_9.block<3,1>(0,3).cross( M10_9.block<3,1>(0,2) );

        Jfull_10.block<3,1>(3,0) = M10_0.block<3,1>(0,0);
        Jfull_10.block<3,1>(3,1) = M10_1.block<3,1>(0,1);
        Jfull_10.block<3,1>(3,2) = M10_2.block<3,1>(0,2);
        Jfull_10.block<3,1>(3,3) = M10_3.block<3,1>(0,1);
        Jfull_10.block<3,1>(3,4) = M10_4.block<3,1>(0,2);
        Jfull_10.block<3,1>(3,5) = M10_5.block<3,1>(0,0);
        Jfull_10.block<3,1>(3,6) = M10_6.block<3,1>(0,1);
        Jfull_10.block<3,1>(3,7) = M10_7.block<3,1>(0,0);
        Jfull_10.block<3,1>(3,8) = M10_8.block<3,1>(0,1);
        Jfull_10.block<3,1>(3,9) = M10_9.block<3,1>(0,2);

        Eigen::Matrix<double, 6, 6> rot = Eigen::Matrix<double, 6, 6>::Zero() ;
        rot.block<3,3>(0,0) = M010.block<3,3>(0,0);
        rot.block<3,3>(3,3) = M010.block<3,3>(0,0);
        Eigen::Matrix<double, 6, 10> Jfull_0= rot*Jfull_10;

        return Jfull_0;
    }

    Vector10d secondTaskEta(const Vector10d &utheta) {
        // to minimize the square of angs of waist.
        // theta(1~3)=0  ---> utheta(1~3) = 0, 0.797473, 0
        Vector10d eta = Vector10d::Zero();
        eta.block<3,1>(0, 0) = 2*Eigen::Vector3d(utheta(0), utheta(1)-0.797473, utheta(2));
        return eta;
    }

    bool checkValid(const Vector10d &angs){
        Vector10d v1 =  angs - cstr.block<10,1>(0,0);
        Vector10d v2 = -angs + cstr.block<10,1>(0,1);
        return v1.minCoeff()>0 && v2.minCoeff()>0;
    }

public:
    void setOptimal(bool opt){
        optFlag = opt;
    }

    void gingerFK(const Vector10d &angs, Eigen::Matrix3d &r, Eigen::Vector3d &v){
        setM(angs);
        r = M010.block<3,3>(0,0);
        v = M010.block<3,1>(0,3);
    }

    bool gingerIK(Vector10d &res, const Eigen::Vector3d &pos, const Eigen::Vector3d &aa, const Vector10d &initialGuess = Vector10d::Zero()){
        int iter = MAX_ITER;
        Eigen::Matrix3d TgtR = Eigen::AngleAxisd(aa.norm(), aa/aa.norm()).toRotationMatrix();
        Eigen::Matrix4d TgtT = Eigen::Matrix4d::Identity();
        TgtT.block<3,3>(0,0) = TgtR;
        TgtT.block<3,1>(0,3) = pos;

        Vector10d curTheta = initialGuess;
    
        Eigen::Matrix3d curR;
        Eigen::Vector3d curV;
        gingerFK(curTheta, curR, curV);

        Eigen::Vector3d v = pos-curV;
        Eigen::Matrix3d deltaR = curR.transpose() * TgtR;
        Eigen::AngleAxisd deltaAA(deltaR);
        Eigen::Vector3d w = deltaAA.axis() * deltaAA.angle(); 
        Eigen::Matrix<double, 6, 1> V(0,0,0,0,0,0);
        V.block<3,1>(0,0) = v;
        V.block<3,1>(3,0) = curR * w;

        while ( --iter>0 && (v.norm()>0.001 || w.norm()>0.01)){
            Eigen::Matrix<double, 6, 10> J = gingerJ(curTheta);
            // if ( condNumber(J) > 10 ){
            //     std::cout<<"WARNING: J is ill-condition"<<std::endl;
            // }

            Eigen::Matrix<double,10,6> pinvJ = J.completeOrthogonalDecomposition().pseudoInverse();
            Vector10d eta = curTheta - (cstr.block<10,1>(0,0)+cstr.block<10,1>(0,1))/2.0;
            Vector10d deltaTheta =  pinvJ * V * dt - (Matrix10d::Identity()-pinvJ*J) * eta *0.5;
            // Vector10d deltaTheta =  pinvJ * V * dt;
            curTheta += deltaTheta;

            gingerFK(curTheta, curR, curV);
            v = pos-curV;
            deltaR = curR.transpose() * TgtR;
            deltaAA.fromRotationMatrix(deltaR);
            w = deltaAA.axis() * deltaAA.angle(); 
            V.block<3,1>(0,0) = v;
            V.block<3,1>(3,0) = curR*w;
        }
        res = curTheta;

        Eigen::Matrix4d diff = M010-TgtT;
        // if (iter>0 && diff.maxCoeff()<0.01 && diff.maxCoeff()>-0.01 && checkValid(res)){
        if (iter>0){
            return true;
        } else {
            // std::cout<<"iters: "<<MAX_ITER-iter<<std::endl;
            // std::cout<<M010<<std::endl;
            // std::cout<<TgtT<<std::endl;
            return false;
        }
    }

    bool gingerIK_Constraint(Vector10d &res, const Eigen::Vector3d &pos, const Eigen::Vector3d &aa, const Vector10d &initialGuess = Vector10d::Zero()){
        int iter = MAX_ITER;
        double k = optFlag ? 0.1 : 0.0;

        Eigen::Matrix3d TgtR = Eigen::AngleAxisd(aa.norm(), aa/aa.norm()).toRotationMatrix();
        Eigen::Matrix4d TgtT = Eigen::Matrix4d::Identity();
        TgtT.block<3,3>(0,0) = TgtR;
        TgtT.block<3,1>(0,3) = pos;

        Vector10d curTheta = initialGuess;
        Vector10d curUTheta = theta2UTheta(initialGuess);

        Eigen::Matrix3d curR;
        Eigen::Vector3d curV;
        gingerFK(curTheta, curR, curV);

        Eigen::Vector3d v = pos-curV;
        Eigen::Matrix3d deltaR = curR.transpose() * TgtR;
        Eigen::AngleAxisd deltaAA(deltaR);
        Eigen::Vector3d w = deltaAA.axis() * deltaAA.angle(); 

        Eigen::Matrix<double, 6, 1> V(0,0,0,0,0,0);
        V.block<3,1>(0,0) = v;
        V.block<3,1>(3,0) = curR*w;
        
        while ( --iter>0 && (v.norm()>0.001 || w.norm()>0.01)){
            Eigen::Matrix<double, 6, 10> J = gingerJ(curTheta);
            Matrix10d uJ = uthetaJ(curUTheta); // positive definite diagonal matrix

            Vector10d diag = uJ.diagonal();
            double cond = diag.maxCoeff()/diag.minCoeff();
            if (cond > condThresh){
                curUTheta = curUTheta / 1000.0;
                for (int i=0; i<10; i++){
                    uJ(i,i) = uJ(i,i) > condThresh ? uJ(i,i)/log(cond) : uJ(i,i);
                    uJ(i,i) = uJ(i,i) < 1.0/condThresh ? uJ(i,i)*log(cond) : uJ(i,i);
                }
            }
            J *= uJ;

            Vector10d eta = secondTaskEta(curUTheta);
            Eigen::Matrix<double, 10, 6> pinvJ = J.completeOrthogonalDecomposition().pseudoInverse();

            Vector10d deltaUTheta = pinvJ*V*dt - (Matrix10d::Identity()-pinvJ*J)*eta*k;
            // Vector10d deltaUTheta = pinvJ*V*dt;

            for (auto item:deltaUTheta) {
                if (isnan(item) || isinf(item)){
                    std::cout<<"delteUTheta NAN or INF ERROR"<<std::endl;
                    return false;
                }
            }

            curUTheta += deltaUTheta;
            for (auto &item : curUTheta) {
                item = fabs(item) < 100.0 ? item : mySign(item)*100.0; // truncate UTheta, little influence to the success rate
            }
            curTheta = utheta2Theta(curUTheta);

            gingerFK(curTheta, curR, curV);
            v = pos-curV;
            deltaR = curR.transpose() * TgtR;
            deltaAA.fromRotationMatrix(deltaR);
            w = deltaAA.axis() * deltaAA.angle(); 
            V.block<3,1>(0,0) = v;
            V.block<3,1>(3,0) = curR*w;
        }

        res = curTheta;

        Eigen::Matrix4d diff = M010-TgtT;
        if (iter>0 && diff.maxCoeff()<0.01 && diff.maxCoeff()>-0.01 && checkValid(res)){
            return true;
        } else {
            // std::cout<<"diff: \n"<<diff<<"\n---"<<std::endl;
            return false;
        }
    }

    bool gingerIK_ConstraintWeight(Vector10d &res, const Eigen::Vector3d &pos, const Eigen::Vector3d &aa,
                                    const Vector10d &initialGuess = Vector10d::Zero(), double posW=1.0, double rotW=1.0){
        int iter = MAX_ITER;
        double k = optFlag ? 0.1 : 0.0;
        if( fabs(posW+rotW)<0.0001 ){
            posW = rotW = 1.0;
        }
        posW = posW/(posW + rotW);
        rotW = 1-rotW;

        Eigen::Matrix3d TgtR = Eigen::AngleAxisd(aa.norm(), aa/aa.norm()).toRotationMatrix();
        Eigen::Matrix4d TgtT = Eigen::Matrix4d::Identity();
        TgtT.block<3,3>(0,0) = TgtR;
        TgtT.block<3,1>(0,3) = pos;

        Vector10d curTheta = initialGuess;
        Vector10d curUTheta = theta2UTheta(initialGuess);

        Eigen::Matrix3d curR;
        Eigen::Vector3d curV;
        gingerFK(curTheta, curR, curV);

        Eigen::Vector3d v = pos-curV;
        Eigen::Matrix3d deltaR = curR.transpose() * TgtR;
        Eigen::AngleAxisd deltaAA(deltaR);
        Eigen::Vector3d w = deltaAA.axis() * deltaAA.angle(); 

        Eigen::Matrix<double, 6, 1> V(0,0,0,0,0,0);
        V.block<3,1>(0,0) = posW*v;
        V.block<3,1>(3,0) = rotW*curR*w;
        
        while ( --iter>0 && (v.norm()>0.001 || w.norm()>0.01)){
            Eigen::Matrix<double, 6, 10> J = gingerJ(curTheta);
            Matrix10d uJ = uthetaJ(curUTheta); // positive definite diagonal matrix

            Vector10d diag = uJ.diagonal();
            double cond = diag.maxCoeff()/diag.minCoeff();
            if (cond > condThresh){
                curUTheta = curUTheta / 1000.0;
                for (int i=0; i<10; i++){
                    uJ(i,i) = uJ(i,i) > condThresh ? uJ(i,i)/log(cond) : uJ(i,i);
                    uJ(i,i) = uJ(i,i) < 1.0/condThresh ? uJ(i,i)*log(cond) : uJ(i,i);
                }
            }
            J *= uJ;

            Vector10d eta = secondTaskEta(curUTheta);
            Eigen::Matrix<double, 10, 6> pinvJ = J.completeOrthogonalDecomposition().pseudoInverse();

            Vector10d deltaUTheta = pinvJ*V*dt - (Matrix10d::Identity()-pinvJ*J)*eta*k;

            for (auto item:deltaUTheta) {
                if (isnan(item) || isinf(item)){
                    std::cout<<"delteUTheta NAN or INF ERROR"<<std::endl;
                    return false;
                }
            }

            curUTheta += deltaUTheta;
            for (auto &item : curUTheta) {
                item = fabs(item) < 100.0 ? item : mySign(item)*100.0; // truncate UTheta, little influence to the success rate
            }
            curTheta = utheta2Theta(curUTheta);

            gingerFK(curTheta, curR, curV);
            v = pos-curV;
            deltaR = curR.transpose() * TgtR;
            deltaAA.fromRotationMatrix(deltaR);
            w = deltaAA.axis() * deltaAA.angle(); 
            V.block<3,1>(0,0) = posW*v;
            V.block<3,1>(3,0) = rotW*curR*w;
        }

        res = curTheta;

        return true;
    }
};