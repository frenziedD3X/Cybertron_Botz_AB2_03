function submitForm() {
    const form = document.getElementById('userForm');
    const formData = {
        name: form.querySelector('[name="name"]').value,
        id: form.querySelector('[name="id"]').value,
        prompt: form.querySelector('[name="prompt"]').value
    };

    fetch('/get_data', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(formData)
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            document.getElementById('result').innerText = `Error: ${data.error}`;
        } else {
            document.getElementById('result').innerText = `Result: ${data.result}`;
            document.getElementById('generated_text').innerText = `Generated Text: ${data.generated_text}`;
        }
    })
    .catch(error => console.error('Error:', error));
}