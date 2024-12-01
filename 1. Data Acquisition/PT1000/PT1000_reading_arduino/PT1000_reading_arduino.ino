#include <Adafruit_MAX31865.h>

// Use software SPI: CS, DI, DO, CLK
Adafruit_MAX31865 thermo1 = Adafruit_MAX31865(10, 11, 12, 13);

// The value of the Rref resistor. Use 430.0 for PT100 and 4300.0 for PT1000
#define RREF      4300.0
// The 'nominal' 0-degrees-C resistance of the sensor
// 100.0 for PT100, 1000.0 for PT1000
#define RNOMINAL  1000.0

void setup() {
  Serial.begin(115200);
  Serial.println("Adafruit MAX31865 PT1000 Sensor Test!");

  thermo1.begin(MAX31865_4WIRE);  // set to 2WIRE or 4WIRE as necessary
}

void loop() {
  // Read temperature from the sensor
  float temperature1 = thermo1.temperature(RNOMINAL, RREF);

  // Print the temperature in a structured format
  Serial.print("TEMP:"); 
  Serial.println(temperature1);

  // Check and print any faults
  uint8_t fault1 = thermo1.readFault();
  if (fault1) {
    Serial.print("FAULT:");
    Serial.println(fault1);
  }

  delay(1000);
}
