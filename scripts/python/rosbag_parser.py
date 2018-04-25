import rosbag, rospy
import yaml, csv, subprocess
import argparse, string
import os
import std_msgs, rosgraph_msgs, sensor_msgs, geometry_msgs
import test_package.msg

if __name__=="__main__":
    #parse inputs
    parser = argparse.ArgumentParser(description="convert rosbags to text files")
    parser.add_argument("rosbag_file", type=str, help="rosbag file name")

    #create a directory to store log files
    current_dir = os.getcwd()
    args = parser.parse_args()
    bagfile = args.rosbag_file
    directory = current_dir + "/" + bagfile[0:len(bagfile)-4]
    try:
        os.makedirs(directory)
    except:
        pass

    print("-------------------------------------------------")
    print("Directory created: ")
    print(directory)
    print("-------------------------------------------------")

    #get bag info
    info_dict = yaml.load(rosbag.Bag(bagfile, "r")._get_yaml_info())
    print("-------------------------------------------------")
    print("Bag summary: ")
    for topic in info_dict["topics"]:
        print(topic)
    print("-------------------------------------------------")

    bag = rosbag.Bag(bagfile)
    bag_content = bag.read_messages()

    #get list of all topics in the bag
    topics = []
    for topic, msg, t in bag_content:
        if topic not in topics:
            topics.append(topic)

    #topics = ["/servo_controls"]
    #process each topic
    for topic in topics:
        print("Processing: {}").format(topic)
        filename = directory + "/" + string.replace(topic, "/" , "_") + ".csv"
        with open(filename,"w") as csvfile:
            filewriter = csv.writer(csvfile, delimiter=';')
            firstIteration = True
            for subtopic, msg, t in bag.read_messages(topic):
                msg_string = str(msg)
                msg_list = string.split(msg_string, "\n")
                data_list = []
                for name_value in msg_list:
                    split_name_value = string.split(name_value, ":")
                    #print(split_name_value)
                    for i in range(len(split_name_value)):
                        split_name_value[i] = string.strip(split_name_value[i])
                        #hack
                    if split_name_value[i] == "[]":
                        split_name_value[i] = ""
                    if split_name_value[i] == "-" or split_name_value[i] == "":
                        pass
                    else:
                        data_list.append(split_name_value)
                #print(data_list)
                if firstIteration:
                    headers = ["time_stamp"]
                    for pair in data_list:
                        headers.append(pair[0])
                    filewriter.writerow(headers)
                    firstIteration = False
                values = [str(t)]
                for pair in data_list:
                    values.append(pair[1])
                filewriter.writerow(values)

    bag.close()