#include "ros/ros.h"
#include "sensor_msgs/JointState.h"

#include "gingerFullIK.hpp"

#include <chrono>
#include <random>

Vector10d signVec(1, -1, 1, -1, 1, -1, -1, -1, -1, 1);
Vector10d initialGuess(0.1, 0.1, 0.1, -0.1, 0.1, 0.1, -0.1, 0, 0, 0);

int main(int argc, char **argv){
    ros::init(argc, argv, "gingerIK");
    ros::NodeHandle nh;
    ros::Publisher jsPub = nh.advertise<sensor_msgs::JointState>("/joint_states", 1);
    sensor_msgs::JointState jsMsg;
    jsMsg.name = {"Back_Z", "Back_X", "Back_Y", "Left_Shoulder_X", "Left_Shoulder_Y", "Left_Elbow_Z",
                "Left_Elbow_X", "Left_Wrist_Z", "Left_Wrist_X", "Left_Wrist_Y"};

    gingerFullIK gik = gingerFullIK();

// #define TEST

#ifdef TEST
    Vector10d angs, res;
    Eigen::Matrix3d R;
    Eigen::Vector3d posTgt;
    Eigen::Matrix4d resM = Eigen::Matrix4d::Identity();

    std::random_device r;
    std::default_random_engine gen{r()};

    int s = 0, v = 0, d = 0;
    int TestNum = 5000;
    Vector10d initial = Vector10d::Ones()*0.01;

    for(int i=0; i<TestNum; i++){
        for (int j=0; j<10; j++){
            // std::uniform_real_distribution<double> distr(gik.cstr(j,0)+0.02, gik.cstr(j,1)-0.02);
            std::uniform_real_distribution<double> distr(gik.cstr(j,0)+0.1, gik.cstr(j,1)-0.1);
            angs(j) = distr(gen);
            // initial(j) = mySign(angs(j))*0.1;
        }
        // initial += angs;
        gik.gingerFK(angs, R, posTgt);

        Eigen::AngleAxisd aa(R);
        Eigen::Vector3d rotTgt = aa.axis() * aa.angle();

        gik.setOptimal(false);
        bool ret1 = gik.gingerIK_Constraint(res, posTgt, rotTgt, initial);
        double n1 = res.block<3,1>(0,0).norm();
        gik.setOptimal(true);
        bool ret2 = gik.gingerIK_Constraint(res, posTgt, rotTgt, initial);
        double n2 = res.block<3,1>(0,0).norm();

        if (ret1) {
            s++;
        } 
        if (ret2) {
            v++;
        }
        if (ret1 && ret2 && n2<=n1){
            d++;
        }
    }
    std::cout << "s: "<<100.0*s/TestNum << "%" << std::endl;
    std::cout << "v: "<<100.0*v/TestNum << "%" << std::endl;
    std::cout << "d: "<<100.0*d/TestNum << "%" << std::endl;

#else
    // // success
    // Eigen::Vector3d posTgt(0.1, -0.2, 0.3);
    // Eigen::Vector3d rotTgt(-0.696043, 1.59329, -1.1621);
    // Vector10d initial = Vector10d::Zero();

    // fail #1
    Eigen::Vector3d posTgt(0.427213, 0.45248, 0.254317);
    Eigen::Vector3d rotTgt(-0.696043, 1.59329, -1.1621);
    Vector10d initial(-1.46703, -0.1751, 0.0336801, 1.73022, 0.462543, -0.0724328, 0.201483, 0.445201, 0.158033, 0.3393);

    // // fail #2
    // Eigen::Vector3d posTgt(0.720082, -0.411935, -0.14202);
    // Eigen::Vector3d rotTgt(0.138177,1.78669, 1.73676);
    // Vector10d initial(1.32845, -0.416591, 0.282587, 2.99968, 1.17953, -0.42358, -0.124288, 0.229169, 0.145362, -0.4779);

    // fail #3
    // Eigen::Vector3d posTgt(0.757367, -0.338946, 0.147557);
    // Eigen::Vector3d rotTgt(0.0233093, 2.79378, -0.0777507);
    // Vector10d initial(1.32845, -0.416591, 0.282587, 2.99968, 1.17953, -0.42358, -0.124288, 0.229169, 0.145362, -0.4779);
    
    // // fail #4
    // Eigen::Vector3d posTgt(0.47593, 0.427813, 0.256124);
    // Eigen::Vector3d rotTgt(-0.443013, 1.6231, -1.51685);
    

    initial += Vector10d::Ones()*0.01;

    clock_t start = clock();
    Vector10d res;
    gik.setOptimal(false);
    bool ret = gik.gingerIK_Constraint(res, posTgt, rotTgt);
    std::cout << "time elapsed: " << (double)(clock() - start)/CLOCKS_PER_SEC << "s" << std::endl;
    if (!ret) {
        std::cout<<"failed"<<std::endl;
        return 0;
    }

    std::cout<<"res: "<<res.transpose()<<std::endl;
    Vector10d angs = res.array() * signVec.array();
    // std::cout<<"angs: "<<angs.transpose()<<std::endl;

    jsMsg.position.clear();
    for (int i=0; i<10; i++){
        jsMsg.position.push_back(angs(i));
    }

    for(int i=0; i<5; i++){
        jsPub.publish(jsMsg);
        ros::Rate(5).sleep();
    }
#endif

    return 0;
}