#include <stdio.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <driver/gpio.h>

#define LED_GPIO 4 // GPIO PIN
void app_main(void) { // note to myself: app_main is for running realtime in the reeboting
    
    // SET GPIO
    gpio_reset_pin(LED_GPIO);
    gpio_set_direction(LED_GPIO, GPIO_MODE_OUTPUT);

    // ON
    gpio_set_level(LED_GPIO, 1);
    printf("LED is ON\n");
    vTaskDelay(pdMS_TO_TICKS(1000)); // wait one sec

    // OFF
    gpio_set_level(LED_GPIO, 0);
    printf("LED is OFF\n");

}