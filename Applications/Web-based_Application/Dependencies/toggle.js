function toggleInputField(fieldId) {
  var field = document.getElementById(fieldId);
  if (fieldId === 'inputField') {
      field.style.display = 'block';
      document.getElementById('fileField').style.display = 'none';
  } else if (fieldId === 'fileField') {
      field.style.display = 'block';
      document.getElementById('inputField').style.display = 'none';
  }
}

function toggleInitInput(inputId){
  var init_input = document.getElementById(inputId);
  if (inputId === 'identical_init') {
    init_input.style.display = 'block';
      document.getElementById('customized_init').style.display = 'none';
      document.getElementById('random_init').style.display = 'none';
  } else if (inputId === 'customized_init') {
    init_input.style.display = 'block';
      document.getElementById('identical_init').style.display = 'none';
      document.getElementById('random_init').style.display = 'none';
  } else if (inputId === 'random_init') {
    init_input.style.display = 'block';
      document.getElementById('identical_init').style.display = 'none';
      document.getElementById('customized_init').style.display = 'none';
  }
}