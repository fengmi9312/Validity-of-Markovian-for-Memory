const ctx = document.getElementById('myChart').getContext('2d');

const labels = ['January', 'February', 'March', 'is $\\theta = \\frac{m}{n}$', 'May', 'June'];
const data = [10, 20, 30, 25, 15, 35];

const chart = new Chart(ctx, {
  type: 'line',
  data: {
    labels: labels,
    datasets: [{
      label: 'My Dataset',
      data: data,
      backgroundColor: 'rgba(0, 123, 255, 0.5)',
      borderColor: 'rgba(0, 123, 255, 1)',
      borderWidth: 1
    }]
  },
  options: {
    title: {
      display: true,
      text: 'This is $\\theta = \\frac{m}{n}$',
      fontSize: 16,
      fontColor: '#333',
      fontFamily: "'Times New Roman', serif",
      padding: 20
    }
  }
});
