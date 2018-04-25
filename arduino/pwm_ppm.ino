#include <PinChangeInterrupt.h>

#include <ros.h>
#include <std_msgs/Int16MultiArray.h>

#define CHANNEL_NUMBER 6            //number of chanels in DX6i
#define CHANNEL_DEFAULT_VALUE 1500  //default servo value
#define THROTTLE_DEFAULT_VALUE 1100  //default throttle value
#define FRAME_LENGTH 22000           //PPM frame length in microseconds 
#define PULSE_LENGTH 300             //pulse length
#define onState 1                    //polarity of the pulses: 1 is positive, 0 is negative
#define sigPin 13                    //PPM signal output pin on the arduino

/*
 * Define pins used to provide RC PWM signal to Arduino
 * Pins 8, 9 and 10 are used since they work on both ATMega328 and 
 * ATMega32u4 board. So this code will work on Uno/Mini/Nano/Micro/Leonardo
 * See PinChangeInterrupt documentation for usable pins on other boards
 */
const byte channel_pin[] = {8, 9, 10, 11};
volatile unsigned long rising_start[] = {0, 0, 0, 0};
volatile long channel_length[] = {0, 0, 0, 0};

// store ppm signals
int ppm[CHANNEL_NUMBER];


// ROS environment
ros::NodeHandle nh; 
std_msgs::Int16MultiArray pwm_msg;

void controlCallback( const std_msgs::Int16MultiArray& msg)
{
   //update first 4 channels in ppm
   for (int i = 0; i < 4; ++i)
   {
       ppm[i] = msg.data[i];
   }
}

ros::Publisher pub("chatter", &pwm_msg);
ros::Subscriber<std_msgs::Int16MultiArray> sub("servo_controls", &controlCallback );

void setup() 
{
  pinMode(channel_pin[0], INPUT);
  pinMode(channel_pin[1], INPUT);
  pinMode(channel_pin[2], INPUT);
  pinMode(channel_pin[3], INPUT);
  
  attachPinChangeInterrupt(digitalPinToPinChangeInterrupt(channel_pin[0]), onRising0, CHANGE);
  attachPinChangeInterrupt(digitalPinToPinChangeInterrupt(channel_pin[1]), onRising1, CHANGE);
  attachPinChangeInterrupt(digitalPinToPinChangeInterrupt(channel_pin[2]), onRising2, CHANGE);
  attachPinChangeInterrupt(digitalPinToPinChangeInterrupt(channel_pin[3]), onRising3, CHANGE);

  // PPM generation set up
  for(int i=0; i<CHANNEL_NUMBER; i++)
  {
      ppm[i]= CHANNEL_DEFAULT_VALUE;
  }

  pinMode(sigPin, OUTPUT);
  digitalWrite(sigPin, !onState);  //set the PPM signal pin to the default state (off)

  noInterrupts();
  TCCR1A = 0;               // set entire TCCR1 register to 0
  TCCR1B = 0;
  OCR1A = 100;              // compare match register, change this
  TCCR1B |= (1 << WGM12);   // turn on CTC mode
  TCCR1B |= (1 << CS11);    // 8 prescaler: 0,5 microseconds at 16mhz
  TIMSK1 |= (1 << OCIE1A);  // enable timer compare interrupt
  interrupts();

  //set up ROS
  nh.getHardware()->setBaud(250000);
  nh.initNode();
  
  //pwm_msg.layout.dim_length = 1;
  pwm_msg.layout.dim = (std_msgs::MultiArrayDimension *)
  malloc(sizeof(std_msgs::MultiArrayDimension) * 1);
  pwm_msg.layout.dim[0].label = "lox";
  pwm_msg.layout.dim[0].size = 4;
  pwm_msg.layout.dim[0].stride = 4;
  pwm_msg.layout.data_offset = 0;
  pwm_msg.data = (int *)malloc(sizeof(int)*4);
  pwm_msg.data_length = 4;

  //initialize pwm_msg channel
  for (int i = 0; i < 4; ++i)
  {
    pwm_msg.data[i] = CHANNEL_DEFAULT_VALUE;
  }

  nh.advertise(pub);
  nh.subscribe(sub);
}


void processPin(byte pin) {
uint8_t trigger = getPinChangeInterruptTrigger(digitalPinToPCINT(channel_pin[pin]));
if(trigger == RISING) {
    rising_start[pin] = micros();
  } else if(trigger == FALLING) {
    channel_length[pin] = micros() - rising_start[pin];
  }
}

void onRising0(void) {
processPin(0);
}
void onRising1(void) {
processPin(1);
}
void onRising2(void) {
processPin(2);
}
void onRising3(void) {
  processPin(3);
}

byte counter = 0;

void loop() 
{
  //noInterrupts();
  for (int i = 0; i < 4; ++i)
  {
    pwm_msg.data[i] = channel_length[i];
  }
  //interrupts();
  if (counter > 1)
  {
      pub.publish( &pwm_msg );
      counter = 0;
  }
  ++counter;
  
  nh.spinOnce();
  delay(10);
}

ISR(TIMER1_COMPA_vect)
{  
  static boolean state = true;
  TCNT1 = 0;
  if (state)
  {  //start pulse
    digitalWrite(sigPin, onState);
    OCR1A = PULSE_LENGTH * 2;
    state = false;
  }
  else
  {  //end pulse and calculate when to start the next pulse
    static byte cur_chan_numb;
    static unsigned int calc_rest;
    digitalWrite(sigPin, !onState);
    state = true;
    if(cur_chan_numb >= CHANNEL_NUMBER)
    {
      cur_chan_numb = 0;
      calc_rest = calc_rest + PULSE_LENGTH; 
      OCR1A = (FRAME_LENGTH - calc_rest) * 2;
      calc_rest = 0;
    }
    else
    {
      OCR1A = (ppm[cur_chan_numb] - PULSE_LENGTH) * 2;
      calc_rest = calc_rest + ppm[cur_chan_numb];
      cur_chan_numb++;
    }     
  }
}
